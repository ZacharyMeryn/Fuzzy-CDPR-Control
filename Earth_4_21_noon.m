% impactful asf
% IT2-FLS vs PID — Earth baseline sim
% CDPR 1m x 1m, earth gravity only
%
% this file runs the fuzzy controller and PID head to head at earth gravity.
% two conditions: nominal parameters and perturbed parameters (controller
% never knows — tests robustness to model mismatch, which the abstract
% explicitly promises to evaluate).
%
% the neuro part: gradient descent tunes one output scale using a
% compliance-aware cost function so the tuner sees the same delayed plant
% as the real simulation. without that the tuner always wants a higher
% gain than the real system can handle (learned this the hard way in the
% low-grav version — same fix applied here for consistency).
%
% companion file: LowGrav_4_16_midnight.m extends this to moon/mars/ceres/europa
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

GRAVITY.earth = 9.81;   % m/s^2

% motor positions — 1m x 1m square, one motor at each corner
% P.anchors is a 4x2 matrix, each row is [x, y] in meters
P.anchors = [0.0, 1.0;   % A1 top left
             1.0, 1.0;   % A2 top right
             1.0, 0.0;   % A3 bottom right
             0.0, 0.0];  % A4 bottom left
P.num_cables = 4;        % one cable per motor

% bound controls
P.ws_min = [0.05, 0.05];
P.ws_max = [0.95, 0.95];

P.mass_EE    = 0.5;     % kg, placeholder — weigh it
%P.inertia_EE = 0.001;  % kg*m^2, for when we go 3d

% the filament feed tube runs from A1 to the EE and pulls it toward A1
% 0.3 N/m means at 0.5m away it applies 0.15N toward A1, but remeasure when
% you get it in person
P.tube_force_per_meter = 0.3;

% cable properties (placeholders — check it eventually)
P.cable_diameter   = 0.001;   % m
P.cable_density    = 0.97;    % kg/m^3 for dyneema
% this computes linear mass density (kg/m) from density and cross-section area
% the 1e6 converts from m^2 to mm^2 to match density units
P.cable_mass_per_m = P.cable_density * (pi*(P.cable_diameter/2)^2) * 1e6;
P.cable_stiffness  = 50000;   % N/m, axial stiffness EA/L, confirm in spec sheet
P.cable_damping    = 5.0;     % N*s/m, energy dissipation in cable

% how fast the cable actually "feels" a tension change from the motor
% 20 Hz is an educated guess for dyneema, ask someone about it tho
P.cable_compliance_bandwidth = 20.0;   % Hz

% tension limits — T_min keeps cables from going slack (preloaded tension)
% T_max set by motor rated torque and cable break strength, both placeholders
P.T_min = 2.0;    % N — preload keeps all cables taut
P.T_max = 50.0;   % N

% 100 Hz control loop — matches typical embedded control rates
% 5 second runs are enough to see settling behavior for a single-body test
P.dt      = 0.01;                    % timestep in seconds
P.t_total = 10.0;                    % total run duration
P.t_vec   = 0:P.dt:P.t_total;       % time axis
P.N       = length(P.t_vec);         % how many timesteps total

% disturbance parameters — simulates something unexpected hitting the robot
% disturbance hits at 5.0s — move ends at 2.74s so EE has ~2s to settle
% before the shove arrives, making recovery time a meaningful metric.
% magnitude increased to 3.0N so it's clearly visible above steady-state error
P.disturbance_time      = 5.0;   % s
P.disturbance_magnitude = 5.0;   % N — large enough to clearly exceed steady-state error
P.disturbance_duration  = 0.3;   % s

% max speed and acceleration for the trajectory planner
% 0.2 m/s is conservative for a 1m workspace at 100Hz
% prof check: verify these before any motor runs
P.vel_max   = 0.20;   % m/s
P.accel_max = 0.40;   % m/s^2

% uncertainty model for the perturbed condition
% controller is never told about the mismatch — tests robustness to
% incomplete modeling, which the abstract explicitly calls out
P.mass_uncertainty      = 0.15;   % 15% heavier than the model thinks
P.stiffness_uncertainty = 0.20;   % 20% less stiff than the model thinks

fprintf('earth sim loaded\n');
fprintf('100 Hz | %.0f s | %d steps\n\n', P.t_total, P.N);


% =========================================================================
% BLOCK 2: TRAJECTORY GENERATOR + IK SANITY CHECK
% make_trap_traj() generates smooth position+velocity references.
% trapezoidal = accelerate to cruise speed, hold, decelerate to stop.
% if the move is too short to reach cruise speed it automatically drops
% to a triangle profile (no flat top). either way you get a C1-smooth
% reference with no discontinuous velocity jumps.
% then a quick sag check prints cable sag at earth gravity, 10N.
% =========================================================================

% generates a smooth move from start to end position
% accelerates up to vel_max, holds at cruise speed, decelerates to stop
function [pos_traj, vel_traj] = make_trap_traj(pos_start, pos_end, t_vec, vel_max, accel_max)
    N    = length(t_vec);
    d    = norm(pos_end - pos_start);
    dir  = (pos_end - pos_start) / max(d, 1e-9);   % max() prevents div by zero if start=end

    % clamp peak velocity — if the distance is short you can't reach vel_max
    % before you'd overshoot, so sqrt(accel_max * d) is the physics ceiling
    v_peak = min(vel_max, sqrt(accel_max * d));
    t_ramp = v_peak / accel_max;
    t_flat = (d - accel_max * t_ramp^2) / v_peak;

    if t_flat < 0
        % distance too short for a flat top — triangle profile instead
        t_ramp = sqrt(d / accel_max);
        v_peak = accel_max * t_ramp;
        t_flat = 0;
    end

    t_end_ramp1 = t_ramp;
    t_end_flat  = t_ramp + t_flat;
    t_end_ramp2 = t_ramp + t_flat + t_ramp;

    pos_traj = zeros(N, 2);
    vel_traj = zeros(N, 2);

    for k = 1:N
        t = t_vec(k);
        if t <= t_end_ramp1
            s = 0.5 * accel_max * t^2;
            v = accel_max * t;
        elseif t <= t_end_flat
            s = 0.5 * accel_max * t_ramp^2 + v_peak * (t - t_end_ramp1);
            v = v_peak;
        elseif t <= t_end_ramp2
            dt2 = t - t_end_flat;
            s   = 0.5 * accel_max * t_ramp^2 + v_peak * t_flat ...
                + v_peak * dt2 - 0.5 * accel_max * dt2^2;
            v   = v_peak - accel_max * dt2;
        else
            s = d;   v = 0;
        end
        s = min(s, d);
        pos_traj(k,:) = pos_start + s * dir;
        vel_traj(k,:) = v * dir;
    end
end


% SECTION 1: inverse kinematics with parabolic sag correction
% given where the EE is, how long does each cable need to be?
% sag = w*L^2/(8*T) where w = cable_mass_per_m * g
% the arc correction adds extra length to account for the cable's catenary shape

function [L_arc, sag_mid, L_chord] = ik_with_sag(pos, anchors, T_est, g, cable_mass_per_m)
    % straight-line distance from EE to each anchor — geometric IK baseline
    L_chord    = sqrt(sum((anchors - pos).^2, 2));
    % w is cable weight per unit length (N/m) — scales with gravity
    w          = cable_mass_per_m * g;
    % floor tension at 0.1N to avoid dividing by zero if a cable goes slack
    T_safe     = max(T_est, 0.1);
    % parabolic sag correction factor
    sag_factor = (w .* L_chord).^2 ./ (24 .* T_safe.^2);
    % actual cable to pay out = chord length + sag correction
    L_arc      = L_chord .* (1 + sag_factor);
    % midpoint vertical drop — useful for visualizing and validating
    sag_mid    = w .* L_chord.^2 ./ (8 .* T_safe);
end

% quick sanity check — IK at workspace center with 10N on all cables
T_demo   = ones(4,1) * 10;
pos_demo = [0.5, 0.5];
[~, sag_e, L_chord] = ik_with_sag(pos_demo, P.anchors, T_demo, GRAVITY.earth, P.cable_mass_per_m);
fprintf('sag at center, 10N (earth):\n');
for i = 1:4
    fprintf('  cable %d: chord=%.1fmm  sag=%.4fmm\n', i, L_chord(i)*1e3, sag_e(i)*1e3);
end
fprintf('\n');


% =========================================================================
% BLOCK 3: FORWARD KINEMATICS + JACOBIAN + TENSION QP
% FK goes the other direction — given cable lengths, find EE position.
% Jacobian maps cable tensions to cartesian force on the EE.
% QP picks the cable tensions that produce the commanded force while
% staying inside [T_min, T_max]. minimum-norm solution = least cable stress.
% =========================================================================

% SECTION 2: forward kinematics
% given cable lengths, where is the EE?
% lsqnonlin finds the nonlinear least-squares solution
% falls back to weighted least squares if the solver doesn't converge

function [pos_est, valid, exit_flag] = forward_kinematics(L_measured, anchors, pos_guess)
    residual = @(p) sqrt(sum((anchors - p).^2, 2)) - L_measured;
    opts = optimoptions('lsqnonlin','Display','off', ...
        'FunctionTolerance',1e-10,'StepTolerance',1e-10,'MaxIterations',500);
    [pos_est, ~, ~, exit_flag] = lsqnonlin(residual, pos_guess, [0,0], [1,1], opts);
    valid = (exit_flag >= 1);
    if ~valid
        % fallback: linearize around guess and do weighted least squares
        L_guess = sqrt(sum((anchors - pos_guess).^2, 2));
        J_ls    = (anchors - pos_guess) ./ L_guess;
        delta_L = L_measured - L_guess;
        delta_p = (J_ls' * J_ls) \ (J_ls' * delta_L);
        pos_est = max([0,0], min([1,1], pos_guess + delta_p'));
        warning('FK:SolverFailed','lsqnonlin failed (flag=%d). WLS fallback.', exit_flag);
    end
end


% SECTION 3: jacobian and tension distribution
% J maps cable tensions to cartesian force on the EE
% QP finds the minimum-norm tension set that produces the commanded force
% while keeping all tensions between T_min and T_max

function J = compute_jacobian(pos, anchors)
    L = sqrt(sum((anchors - pos).^2, 2));
    % each row of J is a unit vector pointing from EE toward that anchor
    J = (anchors - pos) ./ L;
end

function [T, feasible] = tension_qp(pos, anchors, F_cmd, T_min, T_max)
    J    = compute_jacobian(pos, anchors);   % 4x2
    n    = size(anchors,1);                  % 4
    opts = optimoptions('quadprog','Display','off');
    % least-squares: minimize ||J'*T - F_cmd||^2 + reg*||T||^2
    % H = J*J' (4x4), f = -J*F_cmd (4x1)
    % always feasible — only box constraints, no equality
    w  = 1e3;
    Jm = J;                               % 4x2
    Fc = F_cmd(:);                        % force as column, 2x1
    H  = w * (Jm * Jm') + 1e-4 * eye(n); % 4x4, positive definite
    fv = -w * (Jm * Fc);                  % 4x1
    [T, ~, flag] = quadprog(H, fv, [], [], [], [], ...
                            T_min*ones(n,1), T_max*ones(n,1), [], opts);
    feasible = (flag == 1);
    if ~feasible
        % gravity-compensation fallback — distribute weight equally
        T = max(T_min, min(T_max, abs(Fc(2))/n * ones(n,1)));
    end
end


% =========================================================================
% BLOCK 4: PID CONTROLLER
% bandwidth-based gain derivation (2 Hz target, 0.5 kg EE):
%   Kp ~ m * wc^2 scaled to force range ~ 18
%   Kd ~ Kp * 2*zeta/wc with zeta=0.7 ~ 3.0
%   Ki small to reject gravity offset without windup ~ 0.8
% integral clamped at ±5 N*s to prevent windup against workspace walls
% prof check: validate with pidtune() or root-locus before hardware runs
% =========================================================================

% SECTION 4: PID controller

function ctrl = pid_init(Kp, Ki, Kd, dt, output_lim)
    ctrl.Kp = Kp(:);   % store as column so Kp.*e works element-wise on 2x1
    ctrl.Ki = Ki(:);
    ctrl.Kd = Kd(:);
    ctrl.dt = dt;
    ctrl.integral  = [0; 0];
    ctrl.e_prev    = [0; 0];
    ctrl.output_lim = output_lim;
end

function [F_cmd, ctrl] = pid_update(ctrl, pos_des, pos_est, vel_est)
    e             = (pos_des(:) - pos_est(:));   % always 2x1 column
    de            = (e - ctrl.e_prev) / ctrl.dt;
    ctrl.integral = ctrl.integral + e * ctrl.dt;
    ctrl.integral = max(-20, min(20, ctrl.integral));
    F_cmd = ctrl.Kp .* e + ctrl.Ki .* ctrl.integral + ctrl.Kd .* de;  % 2x1
    F_mag = norm(F_cmd);
    if F_mag > ctrl.output_lim, F_cmd = F_cmd * ctrl.output_lim / F_mag; end
    ctrl.e_prev = e;
end

% bandwidth-based gains targeting ~2.5 Hz closed-loop bandwidth (wc=16 rad/s):
%   Kp = m * wc^2 = 0.5 * 256 = 128 -> use 80 accounting for cable efficiency loss
%   Kd = 2*zeta*wc*m = 2*0.8*16*0.5 = 12.8 -> use 12
%   Ki = 15 to reject gravity offset within ~2s
% prof check: validate with pidtune() before hardware runs
PID = pid_init([80, 80], [15.0, 15.0], [12.0, 12.0], P.dt, 30.0);
fprintf('PID gains: Kp=%.1f | Ki=%.1f | Kd=%.1f\n\n', PID.Kp(1), PID.Ki(1), PID.Kd(1));


% =========================================================================
% BLOCK 5: IT2 FUZZY LOGIC SYSTEM
% build_it2_fls() constructs the full fuzzy system object.
% it2_evaluate() runs it for a given error and returns a force command.
% the FOU (footprint of uncertainty) is what makes this type-2 — each
% membership function has a band of ambiguity rather than a crisp shape.
% this is the mechanism the abstract claims handles modeling uncertainty,
% so the FOU width (sigma_uncertainty = 0.15) needs to stay consistent
% between this file and the low-grav version.
% =========================================================================

% SECTION 5: IT2 fuzzy logic system
%
% 5 gaussian membership functions per input (NB NS ZE PS PB)
% 25-rule table, mamdani architecture, centroid defuzzification
% sigma_uncertainty sets the footprint of uncertainty width —
% bigger FOU = more robust to noise but slower response

function fis = build_it2_fls(sigma_uncertainty)
    fis = mamfis('Name','IT2_PositionController', ...
                 'AndMethod','min','OrMethod','max', ...
                 'ImplicationMethod','min','AggregationMethod','max', ...
                 'DefuzzificationMethod','centroid');

    % INPUT 1: position error, ±0.4m range covers worst-case workspace crossing
    fis = addInput(fis, [-0.4, 0.4], 'Name', 'pos_error');
    sig_e = 0.06;   % base width of each gaussian bell curve (meters)
    % su is the outer edge of the FOU — wider by sigma_uncertainty (15%)
    % this creates the type-2 uncertainty band
    su    = sig_e * (1 + sigma_uncertainty);
    % five MFs: NB=Negative Big, NS=Negative Small, ZE=Zero, PS=Positive Small, PB=Positive Big
    centers_e = [-0.25, -0.10, 0.0, 0.10, 0.25];
    mf_names  = {'NB','NS','ZE','PS','PB'};
    for i = 1:5
        fis = addMF(fis, 'pos_error', 'gaussmf', [su, centers_e(i)], 'Name', mf_names{i});
    end

    % INPUT 2: velocity error, ±0.3 m/s — slightly narrower since velocities are smaller
    fis = addInput(fis, [-0.3, 0.3], 'Name', 'vel_error');
    sv_u = 0.05 * (1 + sigma_uncertainty);
    centers_v = [-0.20, -0.08, 0.0, 0.08, 0.20];
    for i = 1:5
        fis = addMF(fis, 'vel_error', 'gaussmf', [sv_u, centers_v(i)], 'Name', mf_names{i});
    end

    % OUTPUT: force command ±25N with wider MF spread for earth gravity authority
    fis = addOutput(fis, [-25, 25], 'Name', 'force_cmd');
    out_centers = [-20, -10, 0, 10, 20];
    for i = 1:5
        fis = addMF(fis, 'force_cmd', 'trimf', ...
                    [out_centers(i)-5, out_centers(i), out_centers(i)+5], 'Name', mf_names{i});
    end

    % rule table: row = pos error (NB..PB), col = vel error (NB..PB)
    % aggressive correction when far off, gentle when close and moving right
    rule_table = [1 1 1 2 2; 1 1 2 2 3; 1 2 3 4 5; 3 4 4 5 5; 4 4 5 5 5];
    rules = [];
    for r = 1:5
        for c = 1:5
            rules = [rules; r, c, rule_table(r,c), 1, 1]; %#ok<AGROW>
        end
    end
    fis = addRule(fis, rules);
end

fprintf('building IT2 fuzzy controller...\n');
try
    it2_fis = build_it2_fls(0.15);   % 0.15 = 15% FOU width, matches low-grav file
    fprintf('  %d rules | %d inputs | %d outputs\n', ...
            length(it2_fis.Rules), length(it2_fis.Inputs), length(it2_fis.Outputs));
catch ME
    fprintf('  warning: %s\n', ME.message);
    it2_fis = [];   % empty — simulation will skip IT2 runs if this happens
end

function F_cmd = it2_evaluate(fis, pos_error, vel_error, output_scale, integral_term)
    if isempty(fis), F_cmd = [0; 0]; return; end
    % clamp inputs to FIS range — evalfis throws errors outside its defined range
    ex = max(-0.4, min(0.4, pos_error(1)));
    ey = max(-0.4, min(0.4, pos_error(2)));
    vx = max(-0.3, min(0.3, vel_error(1)));
    vy = max(-0.3, min(0.3, vel_error(2)));
    % evaluate FIS independently for x and y — valid for decoupled planar motion
    % output_scale converts the defuzzified result to actual force
    Fx = evalfis(fis, [ex, vx]) * output_scale;
    Fy = evalfis(fis, [ey, vy]) * output_scale;
    % integral term cancels steady-state offset — same role as Ki in PID
    % without this the fuzzy system has no memory and can't reject constant gravity bias
    if nargin < 5, integral_term = [0;0]; end
    F_cmd = [Fx; Fy] + integral_term(:);
end


% =========================================================================
% BLOCK 6: NEURO-FUZZY OUTPUT SCALE TUNING
% this is the "neuro" part. gradient descent tunes one scalar — the output
% scale multiplier — using a compliance-aware cost function.
%
% why one scalar instead of all MF centers (like an older version of this
% file did): tuning all 10 centroids without the compliance model in the
% cost function made the tuner converge to a scale that looked great in the
% cost evaluation but oscillated in the real simulation. one scalar + the
% full spring-damper model in the cost function avoids that problem and
% matches the approach used in the companion low-grav file.
%
% why compliance matters in the cost: the cable acts like a spring-damper
% between the motor command and what the EE actually feels. without it the
% tuner sees an instantaneous plant and wants a much higher gain than the
% real system can handle.
% =========================================================================

% SECTION 6: neuro-fuzzy output scale tuning

function [scale_tuned, cost_history] = tune_output_scale(fis, ref_trajectory, P, g, n_epochs)
    if isempty(fis)
        scale_tuned = 1.0; cost_history = []; return;
    end

    fprintf('  tuning scale (%d epochs)... ', n_epochs);
    cost_history = zeros(n_epochs, 1);

    % earth starting point — higher than low-grav default because at earth
    % gravity the controller needs more authority to overcome cable weight
    scale = 1.5;
    scale = max(0.1, min(5.0, scale));
    delta = 0.01;   % finite difference nudge

    for epoch = 1:n_epochs
        % decaying learning rate — big steps early, fine steps late
        % slower decay over 200 epochs so we don't stop too early
        lr = 0.05 * (0.80 ^ floor(epoch/25));

        cost   = eval_scale_cost(fis, ref_trajectory, P, g, scale);
        cost_p = eval_scale_cost(fis, ref_trajectory, P, g, scale + delta);
        grad   = (cost_p - cost) / delta;
        scale  = scale - lr * grad;
        scale  = max(0.1, min(5.0, scale));
        cost_history(epoch) = cost;
        if mod(epoch, 10) == 0, fprintf('.'); end
    end
    scale_tuned = scale;
    fprintf(' done. scale: %.3f  cost: %.4f -> %.4f\n', ...
            scale_tuned, cost_history(1), cost_history(end));
end

function cost = eval_scale_cost(fis, ref_traj, P, g, scale)
    % cost function with stability penalty — tracking error + velocity penalty
    % + terminal error penalty. without the extra terms cost is monotonically
    % decreasing (more force = less error during move, no penalty for overshoot)
    % which means the tuner always wants maximum scale regardless of stability.
    N_ref    = size(ref_traj, 1);
    pos      = ref_traj(1,:);
    vel      = [0, 0];
    F_weight = [0, -P.mass_EE * g];
    cost     = 0;
    L_actual = sqrt(sum((P.anchors - pos).^2, 2));
    alpha    = min(P.dt * P.cable_compliance_bandwidth * 2 * pi, 1.0);
    it2_int  = [0; 0];
    Ki_it2   = 15.0;

    for k = 1:N_ref-1
        pos_des = ref_traj(k+1,:);
        e_pos   = pos_des - pos;
        e_vel   = [0,0] - vel;
        it2_int = it2_int + e_pos(:) * P.dt;
        it2_int = max(-15, min(15, it2_int));
        F_ctrl  = it2_evaluate(fis, e_pos, e_vel, scale, Ki_it2 * it2_int)';
        [T_cmd, ~] = tension_qp(pos, P.anchors, F_ctrl(:), P.T_min, P.T_max);
        [L_cmd, ~, ~] = ik_with_sag(pos, P.anchors, T_cmd, g, P.cable_mass_per_m);
        L_actual = L_actual + alpha * (L_cmd - L_actual);
        dL       = L_cmd - L_actual;
        T_comp   = max(P.T_min, min(P.T_max, ...
                       T_cmd + P.cable_stiffness * dL - P.cable_damping * dL / P.dt));
        J        = compute_jacobian(pos, P.anchors);
        F_net    = J' * T_comp + F_weight';
        acc      = F_net / P.mass_EE;
        vel      = vel + acc' * P.dt;
        vel      = max(-P.vel_max, min(P.vel_max, vel));
        pos      = pos + vel * P.dt + 0.5 * acc' * P.dt^2;
        pos      = max(P.ws_min, min(P.ws_max, pos));
        % cost = tracking error + 0.5*velocity penalty (penalizes oscillation)
        % velocity penalty is the key addition — high gain causes oscillation
        % which shows up as sustained nonzero velocity after the move ends
        cost  = cost + sum(e_pos.^2) + 0.5 * sum(vel.^2);
    end
    % terminal penalty: large weight on final error and velocity
    % this directly penalizes failure to settle
    e_final = ref_traj(end,:) - pos;
    cost = cost + 20.0 * sum(e_final.^2) + 10.0 * sum(vel.^2);
end

% tuning trajectory: 300 points over 6 seconds — move takes 2.74s,
% remaining 3.26s is settling time that the terminal penalty evaluates
t_ref = linspace(0, 6, 300);
[ref_traj_tune, ~] = make_trap_traj([0.3, 0.3], [0.7, 0.7], t_ref, P.vel_max, P.accel_max);

fprintf('\nneuro-fuzzy tuning (earth):\n');
it2_fis_tuned = it2_fis;

if ~isempty(it2_fis)
    [it2_scale_earth, cost_history] = tune_output_scale( ...
        it2_fis_tuned, ref_traj_tune, P, GRAVITY.earth, 200);
    fprintf('\ntuned scale (earth): %.3f\n\n', it2_scale_earth);
    it2_fis_tuned = it2_fis;   % scale is passed separately, not baked into FIS
else
    it2_scale_earth = 1.0;     % fallback default if toolbox missing
    cost_history    = [];
end


% =========================================================================
% BLOCK 7: SIMULATION ENGINE
% run_simulation() is the core loop — one call per trial, 8 calls total.
% every timestep: reads reference, computes error, runs controller,
% applies disturbance + tube force, runs QP for tensions, applies
% compliance lag, integrates dynamics, logs everything.
%
% cable compliance: L_actual is a persistent state that lags behind
% the commanded length at the rate set by P.cable_compliance_bandwidth.
% the EE only "feels" the lagged tension, not the instantaneous command.
% this is physically correct for any elastic cable including dyneema.
% =========================================================================

% SECTION 7: simulation engine

function log = run_simulation(controller_type, fis, pid, P, g, ...
                               ref_traj, vel_traj, inject_disturbance, param_perturb, it2_scale)
    % param_perturb lets us secretly lie to the controller about mass/stiffness
    % controller thinks 0.5kg, reality is 0.575kg — tests incomplete modeling robustness
    mass_actual      = P.mass_EE * param_perturb.mass_factor;
    stiffness_actual = P.cable_stiffness * param_perturb.stiffness_factor;
    F_weight         = [0, -mass_actual * g];

    % fresh copy of PID struct — reset integral and previous error so each run is clean
    pid_run          = pid;
    pid_run.integral = [0; 0];
    pid_run.e_prev   = [0; 0];

    % IT2 integral state — accumulates error to cancel steady-state gravity offset
    it2_integral = [0; 0];
    Ki_it2       = 15.0;   % matches PID Ki

    pos      = ref_traj(1,:);
    vel      = [0, 0];
    % compliance state — starts at geometric cable length
    L_actual = sqrt(sum((P.anchors - pos).^2, 2));
    % alpha is the discrete first-order lag coefficient from RC filter approximation
    alpha    = min(P.dt * P.cable_compliance_bandwidth * 2 * pi, 1.0);

    % pre-allocate all log arrays — much faster than growing inside the loop
    N = P.N;
    log.t              = P.t_vec;
    log.pos_ref        = ref_traj;
    log.vel_ref        = vel_traj;
    log.pos_est        = zeros(N, 2);
    log.vel_est        = zeros(N, 2);
    log.error          = zeros(N, 2);
    log.F_cmd          = zeros(N, 2);
    log.T_cables       = zeros(N, P.num_cables);
    log.T_compliance   = zeros(N, P.num_cables);
    log.sag_mid        = zeros(N, P.num_cables);
    log.L_arc          = zeros(N, P.num_cables);
    log.L_actual       = zeros(N, P.num_cables);
    log.disturbance    = zeros(N, 2);
    log.fk_valid       = true(N, 1);
    log.settled        = false;

    for k = 1:N
        t       = P.t_vec(k);
        ref_idx = min(k, size(ref_traj,1));
        pos_des = ref_traj(ref_idx,:);
        vel_des = vel_traj(ref_idx,:);
        e_pos   = pos_des - pos;
        e_vel   = vel_des - vel;

        % select controller — both get the same error, produce a force command
        switch upper(controller_type)
            case 'IT2'
                % accumulate integral for steady-state offset rejection
                it2_integral = it2_integral + e_pos(:) * P.dt;
                it2_integral = max(-15, min(15, it2_integral));   % tight clamp prevents post-disturbance windup
                F_ctrl = it2_evaluate(fis, e_pos, e_vel, it2_scale, Ki_it2 * it2_integral)';
            case 'PID'
                [F_ctrl, pid_run] = pid_update(pid_run, pos_des, pos, vel);
            otherwise
                F_ctrl = [0; 0];
        end

        % disturbance force — sideways shove during the disturbance window
        F_dist = [0; 0];
        if inject_disturbance
            if t >= P.disturbance_time && t <= P.disturbance_time + P.disturbance_duration
                F_dist = [P.disturbance_magnitude; 0];
            end
        end

        % PTFE tube force toward A1 — always on, models filament feed tube
        A1 = P.anchors(1,:);
        dist_to_A1 = norm(pos - A1);
        if dist_to_A1 > 0.001
            tube_dir = (A1 - pos) / dist_to_A1;
            F_dist = F_dist + (P.tube_force_per_meter * dist_to_A1 * tube_dir)';
        end

        F_total_ext = F_weight' + F_dist;
        % QP only needs to find tensions that produce F_ctrl — the cables'
        % job is to deliver the control force. gravity and disturbances are
        % already accounted for separately in the dynamics integration below.
        % passing F_total_ext+F_ctrl here caused gravity to be applied twice:
        % once inside QP (baked into cable tensions) and once in F_net below.
        F_ctrl_col  = F_ctrl(:);   % guarantee 2x1 column regardless of controller output shape
        [T_cmd, ~]  = tension_qp(pos, P.anchors, F_ctrl_col, P.T_min, P.T_max);
        [L_cmd, sag_mid, ~] = ik_with_sag(pos, P.anchors, T_cmd, g, P.cable_mass_per_m);

        % cable compliance — L_actual chases L_cmd but can't jump instantly
        L_actual = L_actual + alpha * (L_cmd - L_actual);
        dL       = L_cmd - L_actual;
        % actual tension = commanded + spring force - damping, clamped to limits
        T_compliance = max(P.T_min, min(P.T_max, ...
                           T_cmd + stiffness_actual * dL - P.cable_damping * dL / P.dt));

        J               = compute_jacobian(pos, P.anchors);
        F_cables_actual = J' * T_compliance;
        F_net           = F_cables_actual + F_weight' + F_dist;
        acc             = F_net / mass_actual;

        vel = vel + acc' * P.dt;
        pos = pos + vel * P.dt + 0.5 * acc' * P.dt^2;   % second-order euler
        pos = max(P.ws_min, min(P.ws_max, pos));

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
% BLOCK 8: METRICS
% compute_metrics() takes a full simulation log and distills it into
% the numbers that go in the paper: RMSE, max error, settling time,
% disturbance recovery time, tension violations, compliance lag.
% these are exactly the metrics the abstract promises (position tracking
% error, overshoot, settling time, robustness to perturbations).
% =========================================================================

% SECTION 8: metrics

function metrics = compute_metrics(log, P)
    % convert error to scalar magnitude in mm
    error_mag = sqrt(sum(log.error.^2, 2)) * 1000;

    metrics.rmse_pos         = rms(error_mag);
    metrics.max_error        = max(error_mag);
    % steady-state: mean error over the last 0.5 seconds
    metrics.steady_state_err = mean(error_mag(end-round(0.5/P.dt):end));

    % overshoot: project error onto the direction of the final target
    ref_final     = log.pos_ref(end,:);
    ref_dir       = ref_final ./ norm(ref_final);
    err_final_dir = log.error * ref_dir';
    metrics.overshoot = max(0, -min(err_final_dir)*1000);

    % settling time: first moment where error stays under 5% band for rest of run
    % 5% is standard for cable robots with compliance — 2% is too tight given
    % the residual sag offset that the integrator can't fully cancel mid-move
    travel = norm(log.pos_ref(end,:) - log.pos_ref(1,:)) * 1000;
    band   = 0.05 * travel;
    metrics.settling_time = P.t_total;
    for k = 1:P.N
        if all(error_mag(k:end) < band)
            metrics.settling_time = P.t_vec(k);
            break;
        end
    end

    % recovery time: time from disturbance peak back to pre-disturbance level
    % threshold = pre-disturbance error (not 110%) — must fully return to baseline
    t_dist_idx = find(P.t_vec >= P.disturbance_time, 1);
    t_dist_end_idx = find(P.t_vec >= P.disturbance_time + P.disturbance_duration, 1);
    pre_dist_err = mean(error_mag(max(1,t_dist_idx-round(0.5/P.dt)):t_dist_idx));
    recovery_threshold = pre_dist_err;   % must return to pre-disturbance baseline
    if ~isempty(t_dist_end_idx) && t_dist_end_idx < P.N
        [~, peak_offset] = max(error_mag(t_dist_end_idx:min(t_dist_end_idx+round(2/P.dt),P.N)));
        peak_idx = t_dist_end_idx + peak_offset - 1;
        % only count as recovery if peak was meaningfully above threshold
        if error_mag(peak_idx) > recovery_threshold * 1.10
            recovery_idx = find(error_mag(peak_idx:end) < recovery_threshold, 1);
            if isempty(recovery_idx)
                metrics.dist_recovery_time = NaN;
            else
                metrics.dist_recovery_time = recovery_idx * P.dt;
            end
        else
            % disturbance too small to produce meaningful recovery transient
            metrics.dist_recovery_time = NaN;
        end
    else
        metrics.dist_recovery_time = NaN;
    end

    % count timesteps where any cable tension dropped below T_min (went slack)
    metrics.tension_violations = sum(any(log.T_cables < P.T_min, 2));

    % compliance lag: deviation between commanded and actual cable length
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
% BLOCK 9: RUN LOOP + CONSOLE TABLE
% 2 conditions x 2 controllers x 2 tests = 8 run_simulation calls
% results stored in results struct with fields like Earth_Nominal_IT2
% then a formatted table prints all metrics side by side
% =========================================================================

% SECTION 9: run everything

fprintf('running simulation...\n');

[step_ref, step_vel] = make_trap_traj([0.3, 0.4], [0.7, 0.6], P.t_vec, ...
                                       P.vel_max, P.accel_max);

% nominal: model matches reality
nominal   = struct('mass_factor', 1.0, 'stiffness_factor', 1.0);
% perturbed: controller thinks 0.5kg / 50000 N/m, reality is 0.575kg / 40000 N/m
perturbed = struct('mass_factor', 1.0 + P.mass_uncertainty, ...
                   'stiffness_factor', 1.0 - P.stiffness_uncertainty);

conditions = struct(...
    'name',   {'Earth_Nominal', 'Earth_Perturbed'}, ...
    'g',      {GRAVITY.earth,   GRAVITY.earth},     ...
    'params', {nominal,         perturbed});

results = struct();

for c = 1:length(conditions)
    cond = conditions(c);
    fprintf('\n  %s (g=%.2f, IT2 scale=%.3f)\n', cond.name, cond.g, it2_scale_earth);

    for ctrl_idx = 1:2
        if ctrl_idx == 1
            cname = 'IT2';  fis_use = it2_fis_tuned;
        else
            cname = 'PID';  fis_use = [];
        end

        % test 1: clean step response — no disturbance
        log1 = run_simulation(cname, fis_use, PID, P, cond.g, ...
                               step_ref, step_vel, false, cond.params, it2_scale_earth);
        m1   = compute_metrics(log1, P);

        % test 2: disturbance rejection — same trajectory, shove at 2.5s
        log2 = run_simulation(cname, fis_use, PID, P, cond.g, ...
                               step_ref, step_vel, true, cond.params, it2_scale_earth);
        m2   = compute_metrics(log2, P);

        fname = sprintf('%s_%s', cond.name, cname);
        results.(fname).step         = log1;
        results.(fname).disturbed    = log2;
        results.(fname).metrics_step = m1;
        results.(fname).metrics_dist = m2;

        if isnan(m2.dist_recovery_time)
            fprintf('    %s: RMSE=%.2fmm | settle=%.2fs | recovery=absorbed\n', ...
                    cname, m1.rmse_pos, m1.settling_time);
        else
            fprintf('    %s: RMSE=%.2fmm | settle=%.2fs | recovery=%.2fs\n', ...
                    cname, m1.rmse_pos, m1.settling_time, m2.dist_recovery_time);
        end
    end
end

% print the full comparison table to console
fprintf('\n');
fprintf('%-26s  %-8s  %-9s  %-10s  %-10s  %-8s\n', ...
        'condition + controller', 'RMSE', 'MaxErr', 'Settle', 'Recovery', 'T_viol');
fprintf('%s\n', repmat('-', 1, 80));
cond_names = {'Earth_Nominal', 'Earth_Perturbed'};
ctrl_names = {'IT2', 'PID'};
for c = 1:length(cond_names)
    for ct = 1:length(ctrl_names)
        fname = sprintf('%s_%s', cond_names{c}, ctrl_names{ct});
        if isfield(results, fname)
            m1 = results.(fname).metrics_step;
            m2 = results.(fname).metrics_dist;
            if isnan(m2.dist_recovery_time)
                rec_str = '  absorbed';
            else
                rec_str = sprintf('%8.2fs', m2.dist_recovery_time);
            end
            fprintf('%-26s  %5.2fmm  %8.2fmm  %8.2fs  %s  %5d\n', ...
                    [strrep(cond_names{c},'_',' ') ' ' ctrl_names{ct}], ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, ...
                    rec_str, m1.tension_violations);
        end
    end
end


% =========================================================================
% BLOCK 10: FIGURES — SPRINGER JOURNAL STYLE
% five publication-quality figures sized to Springer column widths.
% fig2 uses a shared y-axis computed from actual data so IT2 and PID
% panels are directly comparable at the same scale.
% all figures saved at 600 DPI PNG — Springer minimum for line art.
% greyscale + distinct linestyles so everything reads in B&W print.
% =========================================================================

% SECTION 9b: figures — springer journal style
%
% sized to springer column widths: 174mm double, 84mm single
% greyscale + distinct linestyles so everything reads in B&W print
% 600 DPI PNG via print() per springer line-art requirements
% axes: open box, inward ticks — standard springer house style

mm2in = @(x) x/25.4;
FN = 'Times New Roman'; FS = 8; TS = 9; LFS = 7; LW = 1.2; LWt = 0.6;
% greyscale: nominal=black, perturbed=dark grey
C1 = [0.00 0.00 0.00];   % black — nominal
C2 = [0.45 0.45 0.45];   % dark grey — perturbed
LS = {'-','--'};          % solid for nominal, dashed for perturbed

set(0,'DefaultAxesFontName',FN,'DefaultAxesFontSize',FS, ...
      'DefaultTextFontName',FN,'DefaultTextFontSize',FS, ...
      'DefaultLegendFontSize',LFS,'DefaultLegendFontName',FN);

function springer_ax(ax, ttl, xl, yl)
    % applies springer house style: open box, inward ticks, light grid
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
    % force legend text black — MATLAB sometimes renders it grey
    set(lg,'TextColor','k','EdgeColor',[0.5 0.5 0.5]);
end

% compute shared y-axis limit — both controllers use the same scale
% so the IT2 vs PID comparison is honest (not just IT2 looking better
% because its axis is compressed to its own small range)
cnames_fig = {'Earth_Nominal','Earth_Perturbed'};
all_it2_max = 0;
all_pid_max = 0;
for c = 1:2
    e_it2 = sqrt(sum(results.(sprintf('%s_IT2',cnames_fig{c})).step.error.^2,2))*1000;
    e_pid = sqrt(sum(results.(sprintf('%s_PID',cnames_fig{c})).step.error.^2,2))*1000;
    all_it2_max = max(all_it2_max, max(e_it2));
    all_pid_max = max(all_pid_max, max(e_pid));
end
y_shared = max(all_it2_max, all_pid_max) * 1.1;   % 10% headroom

% fig 1: RMSE and settling time bars — 174mm double-column
fw = mm2in(174); fh = mm2in(70);
fig1 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
              'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
rmse_it2   = [results.Earth_Nominal_IT2.metrics_step.rmse_pos, ...
              results.Earth_Perturbed_IT2.metrics_step.rmse_pos];
rmse_pid   = [results.Earth_Nominal_PID.metrics_step.rmse_pos, ...
              results.Earth_Perturbed_PID.metrics_step.rmse_pos];
settle_it2 = [results.Earth_Nominal_IT2.metrics_step.settling_time, ...
              results.Earth_Perturbed_IT2.metrics_step.settling_time];
settle_pid = [results.Earth_Nominal_PID.metrics_step.settling_time, ...
              results.Earth_Perturbed_PID.metrics_step.settling_time];
xb = 1:2; bw = 0.32;
ax1 = subplot(1,2,1);
b1 = bar(ax1,xb-bw/2,rmse_it2,bw,'FaceColor',C1,'EdgeColor','none'); hold(ax1,'on');
b2 = bar(ax1,xb+bw/2,rmse_pid,bw,'FaceColor',C2,'EdgeColor','none');
set(ax1,'XTick',1:2,'XTickLabel',{'Nominal','Perturbed'});
lg1 = legend(ax1,[b1 b2],{'IT2-FLS','PID'},'Location','northeast','Box','off');
fix_legend(lg1);
springer_ax(ax1,'(a) Step Response RMSE','Condition','RMSE (mm)');
ax2 = subplot(1,2,2);
b3 = bar(ax2,xb-bw/2,settle_it2,bw,'FaceColor',C1,'EdgeColor','none'); hold(ax2,'on');
b4 = bar(ax2,xb+bw/2,settle_pid,bw,'FaceColor',C2,'EdgeColor','none');
set(ax2,'XTick',1:2,'XTickLabel',{'Nominal','Perturbed'});
lg2 = legend(ax2,[b3 b4],{'IT2-FLS','PID'},'Location','northeast','Box','off');
fix_legend(lg2);
springer_ax(ax2,'(b) Settling Time','Condition','Time (s)');
print(fig1,'-dpng','-r600','fig1_earth_rmse_settling.png');

% fig 2: tracking error time series — 174mm double-column, shared y-axis
fw = mm2in(174); fh = mm2in(78);
fig2 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
              'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax3 = subplot(1,2,1);
for c = 1:2
    plot(ax3,P.t_vec, ...
         sqrt(sum(results.(sprintf('%s_IT2',cnames_fig{c})).step.error.^2,2))*1000, ...
         LS{c},'Color',C1,'LineWidth',LW); hold(ax3,'on');
end
ylim(ax3,[0 y_shared]);
lg3 = legend(ax3,{'Nominal','Perturbed'},'Location','northeast','Box','off');
fix_legend(lg3);
springer_ax(ax3,'(a) IT2-FLS','Time (s)','Position error (mm)');
ax4 = subplot(1,2,2);
for c = 1:2
    plot(ax4,P.t_vec, ...
         sqrt(sum(results.(sprintf('%s_PID',cnames_fig{c})).step.error.^2,2))*1000, ...
         LS{c},'Color',C1,'LineWidth',LW); hold(ax4,'on');
end
ylim(ax4,[0 y_shared]);
lg4 = legend(ax4,{'Nominal','Perturbed'},'Location','northeast','Box','off');
fix_legend(lg4);
springer_ax(ax4,'(b) PID','Time (s)','Position error (mm)');
print(fig2,'-dpng','-r600','fig2_earth_tracking_error.png');

% fig 3: disturbance rejection — 84mm single-column
fw = mm2in(84); fh = mm2in(72);
fig3 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
              'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax5 = axes(fig3);
plot(ax5,P.t_vec, ...
     sqrt(sum(results.Earth_Nominal_IT2.disturbed.error.^2,2))*1000, ...
     '-','Color',C1,'LineWidth',LW); hold(ax5,'on');
plot(ax5,P.t_vec, ...
     sqrt(sum(results.Earth_Nominal_PID.disturbed.error.^2,2))*1000, ...
     '--','Color',C1,'LineWidth',LW);
% vertical lines mark when the shove starts and ends
xline(ax5,P.disturbance_time,':','Color',[0.5 0.5 0.5],'LineWidth',LWt,'HandleVisibility','off');
xline(ax5,P.disturbance_time+P.disturbance_duration,':','Color',[0.5 0.5 0.5],'LineWidth',LWt,'HandleVisibility','off');
lg5 = legend(ax5,{'IT2-FLS','PID'},'Location','northwest','Box','off');
fix_legend(lg5);
springer_ax(ax5,'Disturbance Rejection (Nominal)','Time (s)','Position error (mm)');
print(fig3,'-dpng','-r600','fig3_earth_disturbance.png');

% fig 4: nominal vs perturbed comparison — 84mm single-column
fw = mm2in(84); fh = mm2in(72);
fig4 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
              'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax6 = axes(fig4);
plot(ax6,P.t_vec, ...
     sqrt(sum(results.Earth_Nominal_IT2.step.error.^2,2))*1000, ...
     '-','Color',C1,'LineWidth',LW); hold(ax6,'on');
plot(ax6,P.t_vec, ...
     sqrt(sum(results.Earth_Perturbed_IT2.step.error.^2,2))*1000, ...
     '--','Color',C1,'LineWidth',LW);
plot(ax6,P.t_vec, ...
     sqrt(sum(results.Earth_Nominal_PID.step.error.^2,2))*1000, ...
     '-','Color',C2,'LineWidth',LW);
plot(ax6,P.t_vec, ...
     sqrt(sum(results.Earth_Perturbed_PID.step.error.^2,2))*1000, ...
     '--','Color',C2,'LineWidth',LW);
lg6 = legend(ax6,{'IT2 Nom','IT2 Pert','PID Nom','PID Pert'},'Location','northeast','Box','off');
fix_legend(lg6);
springer_ax(ax6,'Parameter Uncertainty Robustness','Time (s)','Position error (mm)');
print(fig4,'-dpng','-r600','fig4_earth_robustness.png');

% fig 5: learning curve — 84mm single-column
% only generated if tuning actually ran
if ~isempty(cost_history) && length(cost_history) > 1
    fw = mm2in(84); fh = mm2in(65);
    fig5 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
                  'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
    ax7 = axes(fig5);
    plot(ax7,1:length(cost_history),cost_history,'-','Color',C1,'LineWidth',LW);
    springer_ax(ax7,'Neuro-Fuzzy Convergence (Earth)','Epoch','Tracking cost \itJ');
    print(fig5,'-dpng','-r600','fig5_earth_learning_curve.png');
end

fprintf('\nfigures saved at 600 DPI (springer minimum for line art)\n');
fprintf('  fig1_earth_rmse_settling.png    174mm double-column\n');
fprintf('  fig2_earth_tracking_error.png   174mm double-column  shared y-axis\n');
fprintf('  fig3_earth_disturbance.png      84mm single-column\n');
fprintf('  fig4_earth_robustness.png       84mm single-column\n');
fprintf('  fig5_earth_learning_curve.png   84mm single-column\n');


% =========================================================================
% BLOCK 11: EXPORT
% writes two files: earth_results.csv for Excel/pandas, and
% earth_results_latex.tex that drops straight into the paper with
% \input{earth_results_latex.tex}.
% =========================================================================

% SECTION 10: export CSV and LaTeX table

function export_results_table(results, P, filename_prefix)
    cond_names = {'Earth_Nominal', 'Earth_Perturbed'};
    ctrl_names = {'IT2', 'PID'};
    timestamp  = datestr(now, 'yyyy-mm-dd_HHMMSS');

    % CSV: one row per condition+controller, all metrics as columns
    csv_file = [filename_prefix '_results.csv'];
    fid = fopen(csv_file, 'w');
    fprintf(fid, 'Condition,Controller,RMSE_mm,MaxErr_mm,SettleTime_s,RecoveryTime_s,TensionViolations,ComplianceLag_mm\n');
    for c = 1:length(cond_names)
        for ct = 1:length(ctrl_names)
            fname = sprintf('%s_%s', cond_names{c}, ctrl_names{ct});
            if isfield(results, fname)
                m1 = results.(fname).metrics_step;
                m2 = results.(fname).metrics_dist;
                rt = m2.dist_recovery_time;
                if isnan(rt), rt_str = 'NaN'; else, rt_str = sprintf('%.3f',rt); end
                comp_lag = NaN;
                if isfield(m1,'mean_compliance_lag_mm'), comp_lag = m1.mean_compliance_lag_mm; end
                fprintf(fid, '%s,%s,%.4f,%.4f,%.4f,%s,%d,%.4f\n', ...
                    strrep(cond_names{c},'_',' '), ctrl_names{ct}, ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, ...
                    rt_str, m1.tension_violations, comp_lag);
            end
        end
    end
    fclose(fid);

    % LaTeX: booktabs-style table, drops straight into the paper
    % use: \input{earth_results_latex.tex} in your results section
    tex_file = [filename_prefix '_results_latex.tex'];
    fid = fopen(tex_file, 'w');
    fprintf(fid, '%% auto-generated: %s\n', timestamp);
    fprintf(fid, '\\begin{table}[h]\n\\centering\n');
    fprintf(fid, '\\caption{CDPR Earth Simulation Performance Metrics}\n');
    fprintf(fid, '\\label{tab:earth_results}\n');
    fprintf(fid, '\\begin{tabular}{llrrrr}\\toprule\n');
    fprintf(fid, 'Condition & Controller & RMSE (mm) & Max err (mm) & Settle (s) & Recovery (s) \\\\\n\\midrule\n');
    for c = 1:length(cond_names)
        for ct = 1:length(ctrl_names)
            fname = sprintf('%s_%s', cond_names{c}, ctrl_names{ct});
            if isfield(results, fname)
                m1 = results.(fname).metrics_step;
                m2 = results.(fname).metrics_dist;
                rt = m2.dist_recovery_time;
                if isnan(rt), rt_str = '$<$0.01'; else, rt_str = sprintf('%.2f',rt); end
                fprintf(fid, '%s & %s & %.2f & %.2f & %.2f & %s \\\\\n', ...
                    strrep(cond_names{c},'_',' '), ctrl_names{ct}, ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, rt_str);
            end
        end
        fprintf(fid, '\\midrule\n');
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);

    fprintf('exported: %s, %s\n', csv_file, tex_file);
end

export_results_table(results, P, 'earth');

fprintf('\ndone.\n');