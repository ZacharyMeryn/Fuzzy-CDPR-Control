% impactful asf
% IT2-FLS vs PID — reduced gravity sim
% CDPR 1m x 1m, four planetary bodies
%
% this file runs the fuzzy controller and PID head to head across
% moon, mars, ceres, and europa. gravity is literally just a number
% that gets swapped out — no hardware changes needed, the physics
% model handles the rest.
%
% the neuro part: gradient descent tunes one output scale per body
% using a compliance-aware cost function so the tuner sees the same
% delayed plant as the real simulation. without that the tuner always
% wants a higher gain than the real system can handle.
%
% companion file: ImpactfulFuzzballs_Earth.m does the earth baseline
%
% FILL THESE IN BEFORE RUNNING FOR REAL:
%   P.mass_EE          weigh the end effector
%   P.cable_diameter   calipers
%   P.cable_mass_per_m dyneema spec sheet
%   P.cable_stiffness  axial stiffness EA/L from spec sheet
%   P.T_min / P.T_max  motor torque rating and cable break strength
%
% needs: Optimization Toolbox (quadprog, lsqnonlin)
%        Fuzzy Logic Toolbox  (R2019b+ for type-2)

% wipe everything from last run — fresh workspace every time
clear; clc; close all;

GRAVITY.earth  = 9.81;
GRAVITY.moon   = 1.62;   
GRAVITY.mars   = 3.72;   
GRAVITY.ceres  = 0.27;   
GRAVITY.europa = 1.32;   

% motor positions — 1m x 1m square, one motor at each corner
% P.anchors is a 4x2 matrix, each row is [x, y] in meters
P.anchors = [0.0, 1.0;   % A1 top left
             1.0, 1.0;   % A2 top right
             1.0, 0.0;   % A3 bottom right
             0.0, 0.0];  % A4 bottom left
P.num_cables = 4;        % one cable per motor, obviously


% bound controls
P.ws_min = [0.05, 0.05];
P.ws_max = [0.95, 0.95];


P.mass_EE    = 0.5;    % kg
%P.inertia_EE = 0.001;  % kg*m^2, for when we go 3d

% the filament feed tube runs from A1 to the EE and pulls it toward A1
% 0.3 N/m means at 0.5m away it applies 0.15N toward A1, but remeasure when
% you get it in person
P.tube_force_per_meter = 0.3;

% cable properties (placeholders, chceck it eventually)
P.cable_diameter   = 0.001;   % m 
P.cable_density    = 0.97;    % kg/m^3 for dyneema
% this computes linear mass density (kg/m) from density and cross-section area
% the 1e6 converts from m^2 to mm^2 to match density units
P.cable_mass_per_m = P.cable_density * (pi*(P.cable_diameter/2)^2) * 1e6;
P.cable_stiffness  = 50000;   % N/m, axial stiffness EA/L, confirm in spec sheet
P.cable_damping    = 5.0;     % N*s/m, energy dissipation in cable

% how fast the cable actually "feels" a tension change from the motor
% 20 Hz is an educated guess for dyneema, ask someone ab it tho
% think of it like: motor commands new length, cable takes .05s or so to
% get there
P.cable_compliance_bandwidth = 20.0;   % Hz

% tension limits, T_min keeps cables from going slack, its preloaded
% tension
% T_max is set by motor rated torque and cable break strength
% both are placeholders, confirm with hardware tho
P.T_min = 2.0;    % N
P.T_max = 50.0;   % N

% 100 Hz control loop — matches typical embedded control rates
% 10 second runs give enough time to see settling behavior
P.dt      = 0.01;                    % timestep in seconds
P.t_total = 10.0;                    % total run duration
P.t_vec   = 0:P.dt:P.t_total;       % time axis — 1001 points for 10s at 100Hz
P.N       = length(P.t_vec);         % how many timesteps total

% disturbance parameters — simulates something unexpected hitting the robot
% at 2.5 seconds a 1.5N sideways force turns on for 0.3 seconds then stops
% the controller never knows it's coming
P.disturbance_time      = 2.5;   % s, when the shove happens
P.disturbance_magnitude = 1.5;   % N, how hard the shove is
P.disturbance_duration  = 0.3;   % s, how long it lasts

% max speed and acceleration for the trajectory planner
% 0.2 m/s is conservative for a 1m workspace at 100Hz
% prof check: verify these before any motor runs
P.vel_max   = 0.20;   % m/s
P.accel_max = 0.40;   % m/s^2

% just prints to console so you know the file loaded without crashing
fprintf('low-grav sim loaded\n');
fprintf('bodies: moon (%.2f) | mars (%.2f) | ceres (%.2f) | europa (%.2f)\n', ...
        GRAVITY.moon, GRAVITY.mars, GRAVITY.ceres, GRAVITY.europa);
fprintf('100 Hz | %.0f s | %d steps\n\n', P.t_total, P.N);


% =========================================================================
% BLOCK 2 (lines ~90-155): TRAJECTORY GENERATOR + IK SANITY CHECK
% make_trap_traj() generates smooth position+velocity references.
% then a quick sag check prints how much cable 1 droops at each body
% under 10N tension at the center — just a sanity check, not used in sim.
% =========================================================================

% generates a smooth move from start to end position
% accelerates up to vel_max, holds at cruise speed, decelerates to stop
% if the distance is too short to reach vel_max it automatically switches
% to a triangle profile (accelerate then immediately decelerate, no flat top)
function [pos_traj, vel_traj] = make_trap_traj(pos_start, pos_end, t_vec, vel_max, accel_max)
    N    = length(t_vec);                            % how many timesteps to fill
    d    = norm(pos_end - pos_start);                % straight-line distance of the move
    dir  = (pos_end - pos_start) / max(d, 1e-9);    % unit vector pointing toward target
                                                     % max(d, 1e-9) prevents div by zero if start=end

    % clamp peak velocity — if the distance is short, you can't reach vel_max
    % before you'd overshoot, so sqrt(accel_max * d) is the physics ceiling
    v_peak = min(vel_max, sqrt(accel_max * d));

    % how long the acceleration ramp takes to reach v_peak
    t_ramp = v_peak / accel_max;
    % how long the flat-top cruise phase lasts
    % formula: total distance minus the two ramp areas, divided by cruise speed
    t_flat = (d - accel_max * t_ramp^2) / v_peak;

    if t_flat < 0
        % distance too short for a flat top — triangle profile instead
        % recalculate ramp time for pure triangle
        t_ramp = sqrt(d / accel_max);
        v_peak = accel_max * t_ramp;
        t_flat = 0;
    end

    % time boundaries for each phase of the profile
    t_end_ramp1 = t_ramp;                  % end of first ramp
    t_end_flat  = t_ramp + t_flat;         % end of cruise phase
    t_end_ramp2 = t_ramp + t_flat + t_ramp; % end of second ramp (= end of move)

    pos_traj = zeros(N, 2);   % pre-allocate position array (N rows, x and y columns)
    vel_traj = zeros(N, 2);   % pre-allocate velocity array

    for k = 1:N
        t = t_vec(k);   % current time
        if t <= t_end_ramp1
            % first ramp — constant acceleration from rest
            s = 0.5 * accel_max * t^2;   % distance traveled: s = 0.5*a*t^2
            v = accel_max * t;            % velocity: v = a*t
        elseif t <= t_end_flat
            % cruise phase — constant velocity, position increases linearly
            s = 0.5 * accel_max * t_ramp^2 + v_peak * (t - t_end_ramp1);
            v = v_peak;
        elseif t <= t_end_ramp2
            % deceleration ramp — symmetric to the first ramp
            dt2 = t - t_end_flat;   % time elapsed since decel started
            s   = 0.5 * accel_max * t_ramp^2 + v_peak * t_flat ...
                + v_peak * dt2 - 0.5 * accel_max * dt2^2;
            v   = v_peak - accel_max * dt2;
        else
            % move complete — hold at destination
            s = d;   v = 0;
        end
        s = min(s, d);                         % clamp so we never overshoot target
        pos_traj(k,:) = pos_start + s * dir;   % actual position = start + distance along direction
        vel_traj(k,:) = v * dir;               % velocity vector = speed * direction
    end
end


% SECTION 1: inverse kinematics with parabolic sag correction
% given where the EE is, how long does each cable need to be?
% sag = w*L^2/(8*T) where w = cable_mass_per_m * g
% this is why lunar sag is worse despite lower gravity —
% lower g means lower minimum tension, which tanks the denominator

function [L_arc, sag_mid, L_chord] = ik_with_sag(pos, anchors, T_est, g, cable_mass_per_m)
    % straight-line distance from EE to each anchor — geometric IK baseline
    L_chord    = sqrt(sum((anchors - pos).^2, 2));
    % w is cable weight per unit length (N/m) — scales with gravity
    w          = cable_mass_per_m * g;
    % floor tension at 0.1N to avoid dividing by zero if a cable goes slack
    T_safe     = max(T_est, 0.1);
    % parabolic sag correction factor — how much extra length the sag adds
    % derived from arc length of a parabola: L_arc = L*(1 + w^2*L^2/(24*T^2))
    sag_factor = (w .* L_chord).^2 ./ (24 .* T_safe.^2);
    % actual cable to pay out = chord length + sag correction
    L_arc      = L_chord .* (1 + sag_factor);
    % midpoint vertical drop — useful for visualizing and validating
    sag_mid    = w .* L_chord.^2 ./ (8 .* T_safe);
end

% quick sanity check — call IK at workspace center with 10N on all cables
% just prints the sag per body to console so you can see the gravity effect
T_demo = ones(4,1) * 10;   % 10N on all four cables, arbitrary test condition
fprintf('sag at center, 10N:\n');
fprintf('  %-8s  %-8s  %-10s\n', 'body', 'g', 'sag (mm)');
bodies_demo = {'moon','mars','ceres','europa'};
g_demo      = [GRAVITY.moon, GRAVITY.mars, GRAVITY.ceres, GRAVITY.europa];
for b = 1:4
    [~, sag_b, ~] = ik_with_sag([0.5,0.5], P.anchors, T_demo, g_demo(b), P.cable_mass_per_m);
    % *1e3 converts meters to mm for readability
    fprintf('  %-8s  %-8.2f  %-10.4f\n', bodies_demo{b}, g_demo(b), sag_b(1)*1e3);
end
fprintf('\n');


% =========================================================================
% BLOCK 3 (lines ~157-243): FORWARD KINEMATICS + JACOBIAN + TENSION QP
% FK goes the other direction — given cable lengths, find EE position.
% Jacobian maps cable tensions to cartesian force.
% QP picks the cable tensions that produce the commanded force while
% staying inside [T_min, T_max].
% =========================================================================

% SECTION 2: forward kinematics
% given cable lengths, where is the EE?
% lsqnonlin finds the nonlinear least-squares solution
% falls back to weighted least squares if the solver doesn't converge

function [pos_est, valid, exit_flag] = forward_kinematics(L_measured, anchors, pos_guess)
    % residual function: how far off are the cable lengths for a candidate position p
    % when residual = 0 everywhere, p is the correct EE position
    residual = @(p) sqrt(sum((anchors - p).^2, 2)) - L_measured;

    % solver options — tight tolerances, lots of iterations
    opts = optimoptions('lsqnonlin','Display','off', ...
        'FunctionTolerance',1e-10,'StepTolerance',1e-10,'MaxIterations',500);

    % run the solver — [0,0] and [1,1] are the lower/upper bounds (workspace limits)
    [pos_est, ~, ~, exit_flag] = lsqnonlin(residual, pos_guess, [0,0], [1,1], opts);

    % exit_flag >= 1 means it converged properly
    valid = (exit_flag >= 1);

    if ~valid
        % fallback: linearize around the initial guess and do a quick least-squares solve
        % less accurate but always gives a finite answer instead of crashing
        L_guess = sqrt(sum((anchors - pos_guess).^2, 2));
        J_ls    = (anchors - pos_guess) ./ L_guess;   % linearized jacobian
        delta_p = (J_ls' * J_ls) \ (J_ls' * (L_measured - L_guess));
        pos_est = max([0,0], min([1,1], pos_guess + delta_p'));
        warning('FK:SolverFailed','lsqnonlin failed (flag=%d). WLS fallback.', exit_flag);
    end
end


% SECTION 3: jacobian and tension distribution
% J maps cable tensions to cartesian force on the EE
% QP finds the minimum-norm tension set that produces the commanded force
% while keeping all tensions between T_min and T_max

function J = compute_jacobian(pos, anchors)
    % L is the length of each cable at current position
    L = sqrt(sum((anchors - pos).^2, 2));
    % each row of J is a unit vector pointing from EE toward that anchor
    % this tells you which direction each cable is pulling
    J = (anchors - pos) ./ L;
end

function [T, feasible] = tension_qp(pos, anchors, F_ext, T_min, T_max)
    J    = compute_jacobian(pos, anchors);
    n    = size(anchors,1);   % number of cables (4)
    opts = optimoptions('quadprog','Display','off');
    % quadprog minimizes 0.5*T'*T (minimize cable tensions)
    % subject to J'*T = F_ext (force balance must hold)
    % and T_min <= T <= T_max (physical cable limits)
    [T, ~, flag] = quadprog(eye(n), zeros(n,1), [], [], J', F_ext(:), ...
                            T_min*ones(n,1), T_max*ones(n,1), [], opts);
    feasible = (flag == 1);   % flag=1 means QP solved cleanly
    % if QP fails (bad geometry, infeasible), default to midrange tensions
    if ~feasible, T = (T_min + T_max)/2 * ones(n,1); end
end


% =========================================================================
% BLOCK 4 (lines ~199-243): PID CONTROLLER
% three functions: pid_init sets up the struct, pid_update runs one step,
% make_pid_for_gravity builds a gain-scheduled version for each body.
% gains scale with sqrt(g/g_earth) — less gravity = lower cable tensions
% = the plant behaves differently so gains need to adjust.
% =========================================================================

% SECTION 4: PID controller
%
% bandwidth-based gain derivation (2 Hz target, 0.5 kg EE):
%   Kp ~ m * wc^2 scaled to force range ~ 50 base
%   Kd ~ Kp * 2*zeta/wc with zeta=0.7
%   Ki small to reject gravity offset without windup
%
% gain scheduling: lower gravity = less cable restoring force = lower gains
% scaling law is sqrt(g/g_earth) — physically motivated but not formally derived
%
% prof check: validate with pidtune() or root-locus before first hardware run

function ctrl = pid_init(Kp, Ki, Kd, dt, output_lim)
    % stores all PID parameters in a struct so we can pass it around easily
    ctrl.Kp = Kp;  ctrl.Ki = Ki;  ctrl.Kd = Kd;
    ctrl.dt = dt;
    ctrl.integral  = [0; 0];    % accumulated error — starts at zero each run
    ctrl.e_prev    = [0; 0];    % last timestep's error — needed for derivative
    ctrl.output_lim = output_lim;  % maximum force magnitude the controller can command
end

function [F_cmd, ctrl] = pid_update(ctrl, pos_des, pos_est, vel_est)
    % error = where you want to be minus where you are (in x and y)
    e  = (pos_des - pos_est)';
    % derivative of error = how fast error is changing this timestep
    de = (e - ctrl.e_prev) / ctrl.dt;

    % accumulate error over time for the I term
    ctrl.integral = ctrl.integral + e * ctrl.dt;
    % clamp integral to ±5 N*s — prevents windup if EE is stuck against a wall
    ctrl.integral = max(-5, min(5, ctrl.integral));

    % the actual PID equation: P + I + D
    F_cmd = ctrl.Kp * e + ctrl.Ki * ctrl.integral + ctrl.Kd * de;

    % limit force magnitude — can't command more than the cables can physically deliver
    F_mag = norm(F_cmd);
    if F_mag > ctrl.output_lim, F_cmd = F_cmd * ctrl.output_lim / F_mag; end

    % save this error for next timestep's derivative calculation
    ctrl.e_prev = e;
end

function pid = make_pid_for_gravity(g, g_earth, dt)
    % scale factor: moon at 1.62 m/s^2 gets scale = sqrt(1.62/9.81) = 0.41
    % so moon PID gains are about 41% of earth gains
    scale = sqrt(g / g_earth);
    pid = pid_init([50*scale, 50*scale], [2.0*scale, 2.0*scale], ...
                   [8.0*scale, 8.0*scale], dt, 20.0);
end

% build one PID struct per body — these get passed to run_simulation later
PID_earth  = make_pid_for_gravity(GRAVITY.earth,  GRAVITY.earth, P.dt);
PID_moon   = make_pid_for_gravity(GRAVITY.moon,   GRAVITY.earth, P.dt);
PID_mars   = make_pid_for_gravity(GRAVITY.mars,   GRAVITY.earth, P.dt);
PID_ceres  = make_pid_for_gravity(GRAVITY.ceres,  GRAVITY.earth, P.dt);
PID_europa = make_pid_for_gravity(GRAVITY.europa, GRAVITY.earth, P.dt);

fprintf('PID gains (Kp): earth=%.2f | moon=%.2f | mars=%.2f | ceres=%.2f | europa=%.2f\n\n', ...
        PID_earth.Kp(1), PID_moon.Kp(1), PID_mars.Kp(1), PID_ceres.Kp(1), PID_europa.Kp(1));


% =========================================================================
% BLOCK 5 (lines ~246-315): IT2 FUZZY LOGIC SYSTEM
% build_it2_fls() constructs the full fuzzy system object.
% it2_evaluate() runs it for a given error and returns a force command.
% the FOU (footprint of uncertainty) is what makes this type-2 — each
% membership function has a band of ambiguity rather than a crisp shape.
% this is the mechanism that handles unquantifiable noise like lunar dust.
% =========================================================================

% SECTION 5: IT2 fuzzy logic system
%
% 5 gaussian membership functions per input (NB NS ZE PS PB)
% 25-rule table, mamdani architecture, centroid defuzzification
% sigma_uncertainty sets the footprint of uncertainty width —
% bigger FOU = more robust to noise but slower response
% this is the mechanism that makes type-2 better than type-1
% for environments with unquantifiable disturbances like lunar dust

function fis = build_it2_fls(sigma_uncertainty)
    % mamfis creates a Mamdani fuzzy inference system object
    % 'min' for AND means take the smaller of two membership values
    % 'centroid' defuzzification means the output is the center of mass of the fuzzy region
    fis = mamfis('Name','IT2_PositionController', ...
                 'AndMethod','min','OrMethod','max', ...
                 'ImplicationMethod','min','AggregationMethod','max', ...
                 'DefuzzificationMethod','centroid');

    % INPUT 1: position error, ±0.4m range covers worst-case crossing of the workspace
    fis = addInput(fis, [-0.4, 0.4], 'Name', 'pos_error');
    sig_e = 0.06;   % base width of each gaussian bell curve (in meters)
    % su is the "upper sigma" — the outer edge of the footprint of uncertainty
    % making it wider by sigma_uncertainty (15%) creates the type-2 uncertainty band
    su    = sig_e * (1 + sigma_uncertainty);
    % five membership functions spread across the error range
    % NB=Negative Big, NS=Negative Small, ZE=Zero, PS=Positive Small, PB=Positive Big
    centers_e = [-0.25, -0.10, 0.0, 0.10, 0.25];
    mf_names  = {'NB','NS','ZE','PS','PB'};
    for i = 1:5
        % gaussmf takes [sigma, center] — each MF is a gaussian at its center
        fis = addMF(fis, 'pos_error', 'gaussmf', [su, centers_e(i)], 'Name', mf_names{i});
    end

    % INPUT 2: velocity error, ±0.3 m/s — slightly narrower since velocities are smaller
    fis = addInput(fis, [-0.3, 0.3], 'Name', 'vel_error');
    sv_u = 0.05 * (1 + sigma_uncertainty);   % same FOU widening logic
    centers_v = [-0.20, -0.08, 0.0, 0.08, 0.20];
    for i = 1:5
        fis = addMF(fis, 'vel_error', 'gaussmf', [sv_u, centers_v(i)], 'Name', mf_names{i});
    end

    % OUTPUT: force command ±20N — triangle MFs on output are standard for Mamdani
    % trimf takes [left_edge, peak, right_edge]
    fis = addOutput(fis, [-20, 20], 'Name', 'force_cmd');
    out_centers = [-16, -8, 0, 8, 16];
    for i = 1:5
        fis = addMF(fis, 'force_cmd', 'trimf', ...
                    [out_centers(i)-4, out_centers(i), out_centers(i)+4], 'Name', mf_names{i});
    end

    % rule table: row = pos error (NB..PB), col = vel error (NB..PB)
    % entry = which output MF to fire (1=NB force, 3=ZE, 5=PB force)
    % aggressive correction when far off, gentle when close and moving right
    rule_table = [1 1 1 2 2; 1 1 2 2 3; 1 2 3 4 5; 3 4 4 5 5; 4 4 5 5 5];
    rules = [];
    for r = 1:5
        for c = 1:5
            % each rule row: [input1_MF, input2_MF, output_MF, weight, AND/OR]
            rules = [rules; r, c, rule_table(r,c), 1, 1]; %#ok<AGROW>
        end
    end
    fis = addRule(fis, rules);
end

% build the FIS — wrapped in try/catch in case Fuzzy Logic Toolbox isn't installed
fprintf('building IT2 fuzzy controller...\n');
try
    it2_fis = build_it2_fls(0.15);   % 0.15 = 15% FOU width
    fprintf('  %d rules | %d inputs | %d outputs\n', ...
            length(it2_fis.Rules), length(it2_fis.Inputs), length(it2_fis.Outputs));
catch ME
    fprintf('  warning: %s\n', ME.message);
    it2_fis = [];   % empty — simulation will skip IT2 runs if this happens
end

function F_cmd = it2_evaluate(fis, pos_error, vel_error, output_scale)
    if isempty(fis), F_cmd = [0; 0]; return; end   % toolbox not available, output nothing
    % clamp inputs to the FIS range — evalfis throws errors outside its defined range
    ex = max(-0.4, min(0.4, pos_error(1)));
    ey = max(-0.4, min(0.4, pos_error(2)));
    vx = max(-0.3, min(0.3, vel_error(1)));
    vy = max(-0.3, min(0.3, vel_error(2)));
    % evaluate FIS independently for x and y — valid for decoupled planar motion
    % output_scale multiplies the defuzzified result to get the actual force
    Fx = evalfis(fis, [ex, vx]) * output_scale;
    Fy = evalfis(fis, [ey, vy]) * output_scale;
    F_cmd = [Fx; Fy];
end


% =========================================================================
% BLOCK 6 (lines ~318-422): NEURO-FUZZY TUNING
% this is the "neuro" part. gradient descent tunes one scalar per body —
% the output scale multiplier. different gravity environments need different
% scales because cable restoring force scales with g.
% the key fix that made previous versions oscillate: eval_scale_cost now
% includes the full compliance spring-damper. without it the tuner saw an
% instantaneous plant and converged to scales of ~8-10 which caused
% closed-loop oscillation in the actual simulation.
% =========================================================================

% SECTION 6: neuro-fuzzy output scale tuning
%
% this is the "neuro" component — gradient descent on one parameter
% per gravity environment. the parameter is the output scale multiplier
% that converts the defuzzified FIS output to an actual force command.
%
% why this works as neuro-fuzzy:
%   different gravity bodies need different controller gains because
%   the available cable restoring force scales with g. this tunes
%   that gain from trajectory data rather than hand-picking it.
%
% key fix that made it work: the cost function includes the full
% cable compliance model. without that, the tuner sees an instantaneous
% plant and wants a much higher gain than the real system can handle,
% which is why previous versions oscillated at scale ~8-10.

function [scale_tuned, cost_history] = tune_output_scale(fis, ref_trajectory, P, g, n_epochs)
    if isempty(fis)
        % toolbox not available — return a reasonable default and bail
        scale_tuned = 2.0; cost_history = []; return;
    end

    fprintf('  tuning scale (g=%.2f, %d epochs)... ', g, n_epochs);
    cost_history = zeros(n_epochs, 1);   % pre-allocate cost log for learning curve plot

    % gravity-dependent starting point — lower g bodies need higher scale
    % because cables have less restoring force to work with
    % formula: moon at g=1.62 starts at 0.8 + 1.5*(1 - 1.62/9.81) = 2.05
    %          ceres at g=0.27 starts at 0.8 + 1.5*(1 - 0.27/9.81) = 2.26
    g_earth = 9.81;
    scale   = 0.8 + 1.5 * (1 - g/g_earth);
    scale   = max(0.2, min(3.0, scale));   % clamp to safe range
    delta   = 0.01;   % tiny nudge for finite difference gradient estimation

    for epoch = 1:n_epochs
        % decaying learning rate — take big steps early when far from optimum,
        % small steps late for fine convergence. halves roughly every 15 epochs.
        lr = 0.08 * (0.65 ^ floor(epoch/15));

        % evaluate cost at current scale and at scale + delta
        cost   = eval_scale_cost(fis, ref_trajectory, P, g, scale);
        cost_p = eval_scale_cost(fis, ref_trajectory, P, g, scale + delta);
        % finite difference gradient: if nudging right increases cost, go left
        grad   = (cost_p - cost) / delta;
        % gradient step — move scale in the direction that reduces cost
        scale  = scale - lr * grad;
        % hard cap at 3.0 — above this the controller overdrives the compliance
        % lag and starts oscillating
        scale  = max(0.2, min(3.0, scale));
        cost_history(epoch) = cost;
        if mod(epoch, 10) == 0, fprintf('.'); end   % progress dots
    end
    scale_tuned = scale;
    fprintf(' done. scale: %.3f  cost: %.4f -> %.4f\n', ...
            scale_tuned, cost_history(1), cost_history(end));
end

function cost = eval_scale_cost(fis, ref_traj, P, g, scale)
    % compliance-aware cost — runs the full spring-damper model
    % so the tuner finds a scale that's actually stable in simulation
    N_ref    = size(ref_traj, 1);
    pos      = ref_traj(1,:);   % start at beginning of reference trajectory
    vel      = [0, 0];          % EE starts at rest
    F_weight = [0, -P.mass_EE * g];   % gravity pulling EE down
    cost     = 0;
    % initialize compliance state — same as run_simulation
    L_actual = sqrt(sum((P.anchors - pos).^2, 2));
    alpha    = min(P.dt * P.cable_compliance_bandwidth * 2 * pi, 1.0);

    for k = 1:N_ref-1
        pos_des = ref_traj(k+1,:);
        e_pos   = pos_des - pos;
        e_vel   = [0,0] - vel;
        % run the fuzzy controller at the current scale
        F_ctrl  = it2_evaluate(fis, e_pos, e_vel, scale)';
        % combine gravity + controller force, run QP to get cable tensions
        F_ext   = F_weight' + F_ctrl(:);
        [T_cmd, ~] = tension_qp(pos, P.anchors, F_ext, P.T_min, P.T_max);
        % compute commanded cable lengths with sag correction
        [L_cmd, ~, ~] = ik_with_sag(pos, P.anchors, T_cmd, g, P.cable_mass_per_m);
        % cable compliance lag — L_actual chases L_cmd but can't jump instantly
        L_actual = L_actual + alpha * (L_cmd - L_actual);
        dL       = L_cmd - L_actual;   % how much stretch is still remaining
        % tension the EE actually feels = commanded + spring - damper
        T_comp   = max(P.T_min, min(P.T_max, ...
                       T_cmd + P.cable_stiffness * dL - P.cable_damping * dL / P.dt));
        J     = compute_jacobian(pos, P.anchors);
        F_net = J' * T_comp + F_weight';   % total force on EE
        acc   = F_net / P.mass_EE;
        vel   = vel + acc' * P.dt;
        vel   = max(-P.vel_max, min(P.vel_max, vel));   % velocity clamp
        pos   = pos + vel * P.dt + 0.5 * acc' * P.dt^2; % second-order euler
        pos   = max(P.ws_min, min(P.ws_max, pos));       % workspace clamp
        % add squared error to cost — we're minimizing total tracking error
        cost  = cost + sum(e_pos.^2);
    end
end

% 300-point tuning trajectory — diagonal move across workspace
t_ref = linspace(0, 3, 300);
[ref_traj_tune, ~] = make_trap_traj([0.3, 0.3], [0.7, 0.7], t_ref, P.vel_max, P.accel_max);

fprintf('\nneuro-fuzzy tuning (one scale per body):\n');
it2_fis_tuned = it2_fis;   % start from the base FIS, tuning only changes the scale

% body lists for tuning loop
body_names_tune = {'Moon','Mars','Ceres','Europa'};
body_g_tune     = [GRAVITY.moon, GRAVITY.mars, GRAVITY.ceres, GRAVITY.europa];
it2_scales      = zeros(1,4);     % will hold the tuned scale for each body
cost_histories  = cell(1,4);      % will hold cost curve for each body

if ~isempty(it2_fis)
    for b = 1:4
        [it2_scales(b), cost_histories{b}] = tune_output_scale( ...
            it2_fis_tuned, ref_traj_tune, P, body_g_tune(b), 80);
    end
    fprintf('\ntuned scales: moon=%.3f | mars=%.3f | ceres=%.3f | europa=%.3f\n\n', ...
            it2_scales(1), it2_scales(2), it2_scales(3), it2_scales(4));
    cost_history = cost_histories{1};   % moon curve used for learning curve plot
else
    it2_scales   = [1.8, 1.2, 2.2, 1.6];   % fallback defaults if toolbox missing
    cost_history = [];
end


% =========================================================================
% BLOCK 7 (lines ~425-530): SIMULATION ENGINE
% run_simulation() is the core loop — one call per trial, 16 calls total.
% every timestep it: reads reference, computes error, runs controller,
% applies disturbance + tube force, runs QP for tensions, applies
% compliance lag, integrates dynamics, logs everything.
% =========================================================================

% SECTION 7: simulation engine
%
% runs one full trial and logs everything.
%
% cable compliance: L_actual is a persistent state that lags behind
% the commanded length at the rate set by P.cable_compliance_bandwidth.
% the EE only "feels" the lagged tension, not the instantaneous command.
% this is physically correct for any elastic cable including dyneema.
%
% euler integration matches what a real-time control loop would do.

function log = run_simulation(controller_type, fis, pid, P, g, ...
                               ref_traj, vel_traj, inject_disturbance, param_perturb, it2_scale)
    % param_perturb lets us secretly lie to the controller about mass/stiffness
    % used for the perturbed condition — controller thinks 0.5kg, reality is 0.575kg
    mass_actual      = P.mass_EE * param_perturb.mass_factor;
    stiffness_actual = P.cable_stiffness * param_perturb.stiffness_factor;
    % constant downward force from gravity — scales with actual (secret) mass
    F_weight         = [0, -mass_actual * g];

    % fresh copy of PID struct — reset integral and previous error so each run is clean
    pid_run          = pid;
    pid_run.integral = [0; 0];
    pid_run.e_prev   = [0; 0];

    % EE starts at first position in the trajectory, at rest
    pos      = ref_traj(1,:);
    vel      = [0, 0];
    % compliance state — L_actual tracks actual cable length, starts at geometric length
    L_actual = sqrt(sum((P.anchors - pos).^2, 2));
    % alpha is the discrete first-order lag coefficient
    % alpha = dt * bandwidth * 2*pi — derived from RC filter approximation
    alpha    = min(P.dt * P.cable_compliance_bandwidth * 2 * pi, 1.0);

    % pre-allocate all log arrays — much faster than growing them inside the loop
    N = P.N;
    log.t              = P.t_vec;
    log.pos_ref        = ref_traj;
    log.vel_ref        = vel_traj;
    log.pos_est        = zeros(N, 2);           % EE position at each timestep
    log.vel_est        = zeros(N, 2);           % EE velocity at each timestep
    log.error          = zeros(N, 2);           % position error (desired - actual)
    log.F_cmd          = zeros(N, 2);           % force commanded by controller
    log.T_cables       = zeros(N, P.num_cables); % commanded cable tensions
    log.T_compliance   = zeros(N, P.num_cables); % actual tensions after compliance lag
    log.sag_mid        = zeros(N, P.num_cables); % midspan sag per cable
    log.L_arc          = zeros(N, P.num_cables); % commanded cable arc lengths
    log.L_actual       = zeros(N, P.num_cables); % actual cable lengths (compliance state)
    log.disturbance    = zeros(N, 2);           % disturbance force applied
    log.fk_valid       = true(N, 1);            % forward kinematics convergence flag
    log.settled        = false;

    for k = 1:N
        t       = P.t_vec(k);
        ref_idx = min(k, size(ref_traj,1));   % safety clamp on index
        pos_des = ref_traj(ref_idx,:);         % desired position this timestep
        vel_des = vel_traj(ref_idx,:);         % desired velocity (feedforward)
        e_pos   = pos_des - pos;               % position error: where you want to be minus where you are
        e_vel   = vel_des - vel;               % velocity error: desired speed minus actual speed

        % select controller — both get the same error, produce a force command
        switch upper(controller_type)
            case 'IT2'
                % pass the body-specific tuned scale so each gravity gets its own gain
                F_ctrl = it2_evaluate(fis, e_pos, e_vel, it2_scale)';
            case 'PID'
                [F_ctrl, pid_run] = pid_update(pid_run, pos_des, pos, vel);
            otherwise
                F_ctrl = [0; 0];   % unknown controller type — do nothing safely
        end

        % disturbance force — starts at zero, gets a lateral shove during the window
        F_dist = [0; 0];
        if inject_disturbance
            if t >= P.disturbance_time && t <= P.disturbance_time + P.disturbance_duration
                F_dist = [P.disturbance_magnitude; 0];   % sideways push, no vertical
            end
        end

        % PTFE tube force toward A1 anchor — always on, not a disturbance
        % models the filament feed tube that goes along the A1 cable
        % force scales with distance — further from A1 = stronger pull
        A1 = P.anchors(1,:);
        dist_to_A1 = norm(pos - A1);
        if dist_to_A1 > 0.001   % > 1mm so we don't divide by near-zero
            tube_dir = (A1 - pos) / dist_to_A1;   % unit vector pointing toward A1
            F_dist = F_dist + (P.tube_force_per_meter * dist_to_A1 * tube_dir)';
        end

        % combine all external forces and run QP to find cable tensions
        F_total_ext = F_weight' + F_dist;
        % F_combined is total load the cables need to balance (external + control)
        F_combined  = F_total_ext(:) + F_ctrl(:);
        % QP returns the minimum-norm tension set that satisfies force balance
        [T_cmd, ~]  = tension_qp(pos, P.anchors, F_combined, P.T_min, P.T_max);
        % IK gives commanded cable lengths with sag correction
        [L_cmd, sag_mid, ~] = ik_with_sag(pos, P.anchors, T_cmd, g, P.cable_mass_per_m);

        % cable compliance — L_actual chases L_cmd but can't jump instantly
        % alpha controls the chase speed (set by cable_compliance_bandwidth)
        L_actual = L_actual + alpha * (L_cmd - L_actual);
        dL       = L_cmd - L_actual;   % residual stretch still in the cable
        % actual tension = commanded + spring force from stretch - damping from stretch rate
        % clamped to physical tension limits
        T_compliance = max(P.T_min, min(P.T_max, ...
                           T_cmd + stiffness_actual * dL - P.cable_damping * dL / P.dt));

        % compute actual force the EE feels from the compliance tensions (not commanded)
        J               = compute_jacobian(pos, P.anchors);
        F_cables_actual = J' * T_compliance;   % J' maps 4 tensions to 2D force
        F_net           = F_cables_actual + F_weight' + F_dist;   % total force
        acc             = F_net / mass_actual;   % F = ma, solve for a

        % euler integration — update velocity then position
        vel = vel + acc' * P.dt;
        pos = pos + vel * P.dt + 0.5 * acc' * P.dt^2;   % second-order: more accurate
        pos = max(P.ws_min, min(P.ws_max, pos));   % clamp to workspace boundaries

        % log everything for this timestep
        log.pos_est(k,:)      = pos;
        log.vel_est(k,:)      = vel;
        log.error(k,:)        = e_pos;
        log.F_cmd(k,:)        = F_ctrl';
        log.T_cables(k,:)     = T_cmd';
        log.T_compliance(k,:) = T_compliance';
        log.sag_mid(k,:)      = sag_mid';
        log.L_arc(k,:)        = L_cmd';
        log.L_actual(k,:)     = L_actual';
        log.disturbance(k,:)  = F_dist';
    end
end


% =========================================================================
% BLOCK 8 (lines ~533-573): METRICS
% compute_metrics() takes a full simulation log and distills it into
% the numbers that actually go in the paper: RMSE, max error, settling
% time, disturbance recovery time, tension violations, compliance lag.
% =========================================================================

% SECTION 8: metrics

function metrics = compute_metrics(log, P)
    % convert error from x,y components to scalar magnitude, and meters to mm
    error_mag = sqrt(sum(log.error.^2, 2)) * 1000;

    metrics.rmse_pos         = rms(error_mag);                              % headline number
    metrics.max_error        = max(error_mag);                              % worst single timestep
    % steady-state: mean error over the last 0.5 seconds of the run
    metrics.steady_state_err = mean(error_mag(end-round(0.5/P.dt):end));

    % overshoot: project error onto the direction of the final target
    % positive overshoot means EE went past the target at some point
    ref_final     = log.pos_ref(end,:);
    ref_dir       = ref_final ./ norm(ref_final);   % unit vector toward final position
    err_final_dir = log.error * ref_dir';            % error projected onto that direction
    metrics.overshoot = max(0, -min(err_final_dir)*1000);

    % settling time: first moment where error stays under 2% band for rest of run
    travel = norm(log.pos_ref(end,:) - log.pos_ref(1,:)) * 1000;  % total move distance in mm
    band   = 0.02 * travel;                          % 2% of travel distance
    metrics.settling_time = P.t_total;               % default = never settled
    for k = 1:P.N
        if all(error_mag(k:end) < band)   % all remaining errors within band
            metrics.settling_time = P.t_vec(k);
            break;
        end
    end

    % recovery time: how long after the disturbance until error drops below 5mm
    t_dist_idx = find(P.t_vec >= P.disturbance_time, 1);
    if ~isempty(t_dist_idx) && t_dist_idx < P.N
        recovery_idx = find(error_mag(t_dist_idx:end) < 5.0, 1);
        if isempty(recovery_idx)
            metrics.dist_recovery_time = NaN;   % never recovered within run
        else
            metrics.dist_recovery_time = recovery_idx * P.dt;
        end
    else
        metrics.dist_recovery_time = NaN;
    end

    % count timesteps where any cable tension dropped below T_min (went slack)
    % should be zero for a healthy controller
    metrics.tension_violations = sum(any(log.T_cables < P.T_min, 2));

    % compliance lag: mean and max deviation between commanded and actual cable length
    % tells you how much the spring-damper model is affecting the response
    if isfield(log, 'L_actual')
        dL_all = log.L_arc - log.L_actual;
        metrics.mean_compliance_lag_mm = mean(abs(dL_all(:))) * 1000;
        metrics.max_compliance_lag_mm  = max(abs(dL_all(:))) * 1000;
    else
        metrics.mean_compliance_lag_mm = NaN;
        metrics.max_compliance_lag_mm  = NaN;
    end
end


% =========================================================================
% BLOCK 9 (lines ~576-645): RUN LOOP + CONSOLE TABLE
% outer loop: 4 bodies x 2 controllers x 2 tests = 16 run_simulation calls
% each result stored in results struct with field names like Moon_IT2
% then prints a formatted table to console showing all metrics side by side
% =========================================================================

% SECTION 9: run everything
% 4 bodies x 2 controllers x 2 tests = 16 runs total

fprintf('running simulation...\n');

% build the reference trajectory — diagonal move from (0.3,0.4) to (0.7,0.6)
[step_ref, step_vel] = make_trap_traj([0.3, 0.4], [0.7, 0.6], P.t_vec, ...
                                       P.vel_max, P.accel_max);

% nominal: no uncertainty, model matches reality
nominal = struct('mass_factor', 1.0, 'stiffness_factor', 1.0);

% lookup tables mapping body index to name, gravity, and pre-built PID struct
body_names = {'Moon', 'Mars', 'Ceres', 'Europa'};
body_g     = [GRAVITY.moon, GRAVITY.mars, GRAVITY.ceres, GRAVITY.europa];
body_pid   = {PID_moon, PID_mars, PID_ceres, PID_europa};

results = struct();   % will hold all 16 run logs and metrics

for b = 1:length(body_names)
    bname   = body_names{b};
    g       = body_g(b);
    pid_b   = body_pid{b};           % body-specific PID gains
    scale_b = it2_scales(b);         % body-specific IT2 output scale (tuned above)

    fprintf('\n  %s (g=%.2f, IT2 scale=%.3f)\n', bname, g, scale_b);

    for ctrl_idx = 1:2
        if ctrl_idx == 1
            cname = 'IT2';  fis_use = it2_fis_tuned;   % IT2 gets the fuzzy system
        else
            cname = 'PID';  fis_use = [];               % PID doesn't use FIS
        end

        % test 1: clean step response — no disturbance
        log1 = run_simulation(cname, fis_use, pid_b, P, g, ...
                               step_ref, step_vel, false, nominal, scale_b);
        m1   = compute_metrics(log1, P);

        % test 2: disturbance rejection — same trajectory, shove at 2.5s
        log2 = run_simulation(cname, fis_use, pid_b, P, g, ...
                               step_ref, step_vel, true, nominal, scale_b);
        m2   = compute_metrics(log2, P);

        % store results under a field name like 'Moon_IT2' or 'Mars_PID'
        fname = sprintf('%s_%s', bname, cname);
        results.(fname).step         = log1;
        results.(fname).disturbed    = log2;
        results.(fname).metrics_step = m1;
        results.(fname).metrics_dist = m2;
        results.(fname).g            = g;
        results.(fname).it2_scale    = scale_b;

        fprintf('    %s: RMSE=%.2fmm | settle=%.2fs | recovery=%.2fs\n', ...
                cname, m1.rmse_pos, m1.settling_time, m2.dist_recovery_time);
    end
end

% print the full comparison table to console
fprintf('\n');
fprintf('%-22s  %-8s  %-9s  %-10s  %-10s  %-8s\n', ...
        'body + controller', 'RMSE', 'MaxErr', 'Settle', 'Recovery', 'T_viol');
fprintf('%s\n', repmat('-', 1, 75));
ctrl_names = {'IT2','PID'};
for b = 1:length(body_names)
    for ct = 1:2
        fname = sprintf('%s_%s', body_names{b}, ctrl_names{ct});
        if isfield(results, fname)
            m1 = results.(fname).metrics_step;
            m2 = results.(fname).metrics_dist;
            fprintf('%-22s  %5.2fmm  %8.2fmm  %8.2fs  %8.2fs  %5d\n', ...
                    [body_names{b} ' ' ctrl_names{ct}], ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, ...
                    m2.dist_recovery_time, m1.tension_violations);
        end
    end
end


% =========================================================================
% BLOCK 10 (lines ~648-793): FIGURES
% five publication-quality figures sized to Springer column widths.
% key detail: fig2 uses a shared y-axis computed from the actual data
% so the IT2 and PID panels are directly comparable at the same scale.
% all figures saved at 600 DPI PNG — Springer minimum for line art.
% =========================================================================

% SECTION 9b: figures — springer journal style
%
% sized to springer column widths: 174mm double, 84mm single
% greyscale + distinct linestyles so everything reads in B&W print
% 600 DPI PNG via print() per springer line-art requirements
% axes: open box, inward ticks — standard springer house style

mm2in = @(x) x/25.4;   % springer specifies dimensions in mm — convert to inches for MATLAB
% typography settings — all match springer body text specifications
FN = 'Times New Roman'; FS = 8; TS = 9; LFS = 7; LW = 1.2; LWt = 0.6;
% greyscale palette — four shades that are distinguishable in B&W print
C1 = [0.00 0.00 0.00];   % black — Moon
C2 = [0.40 0.40 0.40];   % dark grey — Mars
C3 = [0.70 0.70 0.70];   % light grey — Ceres
C4 = [0.20 0.20 0.20];   % near-black — Europa (different linestyle differentiates it from Moon)
bclr = {C1,C2,C3,C4};   % lookup table, indexed by body
% different linestyle per body — this + color = readable even when photocopied in B&W
LS = {'-','--',':','-.'};

% set as MATLAB global defaults so every figure inherits these automatically
set(0,'DefaultAxesFontName',FN,'DefaultAxesFontSize',FS, ...
      'DefaultTextFontName',FN,'DefaultTextFontSize',FS, ...
      'DefaultLegendFontSize',LFS,'DefaultLegendFontName',FN);

function springer_ax(ax, ttl, xl, yl)
    % applies springer house style to an axis: open box, inward ticks, light grid
    set(ax,'Box','off','TickDir','in','LineWidth',0.5,'Color','white', ...
           'XColor','k','YColor','k','GridColor',[0.85 0.85 0.85], ...
           'GridAlpha',1.0,'FontName','Times New Roman','FontSize',8);
    grid(ax,'on');
    if ~isempty(ttl)
        title(ax,ttl,'FontSize',9,'FontWeight','bold', ...
              'FontName','Times New Roman','Color','k');
    end
    if ~isempty(xl), xlabel(ax,xl,'FontSize',8,'FontName','Times New Roman','Color','k'); end
    if ~isempty(yl), ylabel(ax,yl,'FontSize',8,'FontName','Times New Roman','Color','k'); end
end

function fix_legend(lg)
    % force legend text black — MATLAB sometimes renders it grey by default
    set(lg,'TextColor','k','EdgeColor',[0.5 0.5 0.5]);
end

% compute shared y-axis limit for fig2 — both panels use the same scale
% this makes the IT2 vs PID comparison honest (not just IT2 looking better
% because its axis is compressed to its own small range)
all_it2_max = 0;
all_pid_max = 0;
for b = 1:4
    e_it2 = sqrt(sum(results.(sprintf('%s_IT2',body_names{b})).step.error.^2,2))*1000;
    e_pid = sqrt(sum(results.(sprintf('%s_PID',body_names{b})).step.error.^2,2))*1000;
    all_it2_max = max(all_it2_max, max(e_it2));
    all_pid_max = max(all_pid_max, max(e_pid));
end
y_shared = max(all_it2_max, all_pid_max) * 1.1;   % 10% headroom above maximum

% fig 1: RMSE and settling time bars — 174mm double-column
fw = mm2in(174); fh = mm2in(70);
fig1 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
              'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
% pull the metrics for all 4 bodies and both controllers
rmse_it2 = zeros(1,4); rmse_pid = zeros(1,4);
settle_it2 = zeros(1,4); settle_pid = zeros(1,4);
for b = 1:4
    rmse_it2(b)   = results.(sprintf('%s_IT2',body_names{b})).metrics_step.rmse_pos;
    rmse_pid(b)   = results.(sprintf('%s_PID',body_names{b})).metrics_step.rmse_pos;
    settle_it2(b) = results.(sprintf('%s_IT2',body_names{b})).metrics_step.settling_time;
    settle_pid(b) = results.(sprintf('%s_PID',body_names{b})).metrics_step.settling_time;
end
xb = 1:4; bw = 0.32;   % bar positions and width
ax1 = subplot(1,2,1);
% bars offset by ±bw/2 so IT2 and PID sit side by side for each body
b1 = bar(ax1,xb-bw/2,rmse_it2,bw,'FaceColor',C1,'EdgeColor','none'); hold(ax1,'on');
b2 = bar(ax1,xb+bw/2,rmse_pid,bw,'FaceColor',C2,'EdgeColor','none');
set(ax1,'XTick',1:4,'XTickLabel',body_names);
lg1 = legend(ax1,[b1 b2],{'IT2-FLS','PID'},'Location','northeast','Box','off');
fix_legend(lg1);
springer_ax(ax1,'(a) Step Response RMSE','Planetary body','RMSE (mm)');
ax2 = subplot(1,2,2);
b3 = bar(ax2,xb-bw/2,settle_it2,bw,'FaceColor',C1,'EdgeColor','none'); hold(ax2,'on');
b4 = bar(ax2,xb+bw/2,settle_pid,bw,'FaceColor',C2,'EdgeColor','none');
set(ax2,'XTick',1:4,'XTickLabel',body_names);
lg2 = legend(ax2,[b3 b4],{'IT2-FLS','PID'},'Location','northeast','Box','off');
fix_legend(lg2);
springer_ax(ax2,'(b) Settling Time','Planetary body','Time (s)');
print(fig1,'-dpng','-r600','fig1_rmse_settling.png');

% fig 2: tracking error time series — 174mm double-column, shared y-axis
fw = mm2in(174); fh = mm2in(78);
fig2 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
              'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax3 = subplot(1,2,1);
for b = 1:4
    plot(ax3,P.t_vec,sqrt(sum(results.(sprintf('%s_IT2',body_names{b})).step.error.^2,2))*1000, ...
         LS{b},'Color',bclr{b},'LineWidth',LW); hold(ax3,'on');
end
ylim(ax3,[0 y_shared]);   % shared scale — same limit as PID panel
lg3 = legend(ax3,body_names,'Location','northeast','Box','off');
fix_legend(lg3);
springer_ax(ax3,'(a) IT2-FLS','Time (s)','Position error (mm)');
ax4 = subplot(1,2,2);
for b = 1:4
    plot(ax4,P.t_vec,sqrt(sum(results.(sprintf('%s_PID',body_names{b})).step.error.^2,2))*1000, ...
         LS{b},'Color',bclr{b},'LineWidth',LW); hold(ax4,'on');
end
ylim(ax4,[0 y_shared]);   % same shared scale
lg4 = legend(ax4,body_names,'Location','northeast','Box','off');
fix_legend(lg4);
springer_ax(ax4,'(b) PID','Time (s)','Position error (mm)');
print(fig2,'-dpng','-r600','fig2_tracking_error.png');

% fig 3: disturbance rejection — 84mm single-column
fw = mm2in(84); fh = mm2in(72);
fig3 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
              'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax5 = axes(fig3);
for b = 1:4
    plot(ax5,P.t_vec,sqrt(sum(results.(sprintf('%s_IT2',body_names{b})).disturbed.error.^2,2))*1000, ...
         LS{b},'Color',bclr{b},'LineWidth',LW); hold(ax5,'on');
end
% dotted vertical lines mark when the shove starts and ends
xline(ax5,P.disturbance_time,':','Color',[0.5 0.5 0.5],'LineWidth',LWt,'HandleVisibility','off');
xline(ax5,P.disturbance_time+P.disturbance_duration,':','Color',[0.5 0.5 0.5],'LineWidth',LWt,'HandleVisibility','off');
lg5 = legend(ax5,body_names,'Location','northwest','Box','off');
fix_legend(lg5);
springer_ax(ax5,'IT2-FLS Disturbance Rejection','Time (s)','Position error (mm)');
print(fig3,'-dpng','-r600','fig3_disturbance.png');

% fig 4: cable sag — 84mm single-column, smoothed to remove compliance noise
fw = mm2in(84); fh = mm2in(72);
fig4 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
              'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax6 = axes(fig4);
for b = 1:4
    raw = results.(sprintf('%s_IT2',body_names{b})).step.sag_mid(:,1)*1000;
    % smooth() with window 21 removes the high-frequency compliance spikes
    % the underlying sag trend is what we care about here
    plot(ax6,P.t_vec,smooth(raw,21),LS{b},'Color',bclr{b},'LineWidth',LW); hold(ax6,'on');
end
lg6 = legend(ax6,body_names,'Location','southeast','Box','off');
fix_legend(lg6);
springer_ax(ax6,'Cable 1 Midspan Sag','Time (s)','Sag (mm)');
print(fig4,'-dpng','-r600','fig4_cable_sag.png');

% fig 5: learning curve — 84mm single-column
% only generated if tuning actually ran (cost_history is non-empty)
if ~isempty(cost_history) && length(cost_history) > 1
    fw = mm2in(84); fh = mm2in(65);
    fig5 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
                  'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
    ax7 = axes(fig5);
    plot(ax7,1:length(cost_history),cost_history,'-','Color',C1,'LineWidth',LW);
    springer_ax(ax7,'Neuro-Fuzzy Convergence (Moon)','Epoch','Tracking cost \itJ');
    print(fig5,'-dpng','-r600','fig5_learning_curve.png');
end

fprintf('\nfigures saved at 600 DPI (springer minimum for line art)\n');
fprintf('  fig1_rmse_settling.png   174mm double-column\n');
fprintf('  fig2_tracking_error.png  174mm double-column  shared y-axis\n');
fprintf('  fig3_disturbance.png     84mm single-column\n');
fprintf('  fig4_cable_sag.png       84mm single-column\n');
fprintf('  fig5_learning_curve.png  84mm single-column\n');


% =========================================================================
% BLOCK 11 (lines ~796-857): EXPORT
% writes two files to working directory: lowgrav_results.csv for Excel/
% pandas, and lowgrav_results_latex.tex that drops straight into a paper
% with \\input{lowgrav_results_latex.tex}.
% =========================================================================

% SECTION 10: export CSV and LaTeX table

function export_results_table(results, P, body_names, filename_prefix)
    ctrl_names = {'IT2', 'PID'};
    timestamp  = datestr(now, 'yyyy-mm-dd_HHMMSS');  % for the LaTeX comment header

    % CSV: one row per body+controller combination, all metrics as columns
    csv_file = [filename_prefix '_results.csv'];
    fid = fopen(csv_file, 'w');
    fprintf(fid, 'Body,g_ms2,Controller,RMSE_mm,MaxErr_mm,SettleTime_s,RecoveryTime_s,TensionViolations,ComplianceLag_mm\n');
    for b = 1:length(body_names)
        for ct = 1:length(ctrl_names)
            fname = sprintf('%s_%s', body_names{b}, ctrl_names{ct});
            if isfield(results, fname)
                m1 = results.(fname).metrics_step;
                m2 = results.(fname).metrics_dist;
                g  = results.(fname).g;
                rt = m2.dist_recovery_time;
                % NaN recovery time means it never recovered — write the string 'NaN'
                if isnan(rt), rt_str = 'NaN'; else, rt_str = sprintf('%.3f',rt); end
                comp_lag = NaN;
                if isfield(m1,'mean_compliance_lag_mm'), comp_lag = m1.mean_compliance_lag_mm; end
                fprintf(fid, '%s,%.2f,%s,%.4f,%.4f,%.4f,%s,%d,%.4f\n', ...
                    body_names{b}, g, ctrl_names{ct}, ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, ...
                    rt_str, m1.tension_violations, comp_lag);
            end
        end
    end
    fclose(fid);

    % LaTeX: booktabs-style table, drops straight into a paper
    % use: \input{lowgrav_results_latex.tex} in your results section
    tex_file = [filename_prefix '_results_latex.tex'];
    fid = fopen(tex_file, 'w');
    fprintf(fid, '%% auto-generated: %s\n', timestamp);
    fprintf(fid, '\\begin{table}[h]\n\\centering\n');
    fprintf(fid, '\\caption{CDPR Low-Gravity Simulation Performance Metrics}\n');
    fprintf(fid, '\\label{tab:lowgrav_results}\n');
    fprintf(fid, '\\begin{tabular}{llrrrrr}\n\\toprule\n');
    fprintf(fid, 'Body & Controller & $g$ (m/s\\textsuperscript{2}) & RMSE (mm) & Max err (mm) & Settle (s) & Recovery (s) \\\\\n\\midrule\n');
    for b = 1:length(body_names)
        for ct = 1:length(ctrl_names)
            fname = sprintf('%s_%s', body_names{b}, ctrl_names{ct});
            if isfield(results, fname)
                m1 = results.(fname).metrics_step;
                m2 = results.(fname).metrics_dist;
                g  = results.(fname).g;
                rt = m2.dist_recovery_time;
                if isnan(rt), rt_str = '---'; else, rt_str = sprintf('%.2f',rt); end
                fprintf(fid, '%s & %s & %.2f & %.2f & %.2f & %.2f & %s \\\\\n', ...
                    body_names{b}, ctrl_names{ct}, g, ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, rt_str);
            end
        end
        fprintf(fid, '\\midrule\n');   % horizontal rule between bodies
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);

    fprintf('exported: %s, %s\n', csv_file, tex_file);
end

export_results_table(results, P, body_names, 'lowgrav');

fprintf('\ndone.\n');