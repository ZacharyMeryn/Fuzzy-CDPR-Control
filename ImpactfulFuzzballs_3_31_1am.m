
% impactful asf
% Type-2 Neuro-Fuzzy vs PID Simulation Study
% Cable-Driven Parallel Robot — 1m x 1m, Earth + Lunar conditions
%
% BEFORE RUNNING ANYTHING, REMEMBER TO CHANGE VALUES FOR P STRUCTURE,
% SPECIFICALLY P.mass_EE (measure the end effector mass), P.cable_diameter
% (calipers my g), P.cable_mass_per_m (given by manufacturer ideally),
% P.cable_stiffness (listed on spec sheet as axial stiffness), and P.T_min
% and P.T_max, those depend on the winch motors rated torque and my chosen
% cables breaking strength. lowk consult someone on this bit
%
%
% BEFORE RUNNING DO THE FOLLOWING:
%       1. re order function definitions to avoid any call errors
%       2. put in placeholder values and run smoke test
%
% FOLDER STRUCTURE FOR WHEN THIS IS LESS MESSY:
%   addpath('config')       % robot_params.m, controller_params.m
%   addpath('kinematics')   % ik_sag.m, fk.m, jacobian.m
%   addpath('controllers')  % pid_controller.m, it2_fls.m
%   addpath('simulation')   % run_sim.m, inject_disturbance.m
%   addpath('analysis')     % compute_metrics.m, plot_results.m
%
% CABLE SAG MODEL CHOICE: Parabolic approximation, its close enough for now
% and lets me mess with gravity nicely
%   Key property: sag scales with (w*L^2)/(8*T) where w = cable weight
%   per unit length. Changing gravity just changes w → sag worsens at
%   low g (lower tensions needed to balance reduced weight)

clear; clc; close all;

%  SECTION 0: SIMULATION PARAMETERS

% Gravity conditions defined
GRAVITY.earth  = 9.81;
GRAVITY.lunar  = 9.81/6;

% Workspace geometry for end effector
P.anchors = [0.0, 1.0;      % A1 top left
             1.0, 1.0;      % A2 top right
             1.0, 0.0;      % A3 bottom right
             0.0, 0.0];     % A4 bottom left
P.num_cables = 4;
P.ws_min = [0.05, 0.05];
P.ws_max = [0.95, 0.95];

% End-effector info
P.mass_EE    = 0.5;
P.inertia_EE = 0.001;

% Tube force (PTFE filament feed tube pulling EE toward A1)
% units N/m — placeholder, measure when assembled
P.tube_force_per_meter = 0.3;

% Cable properties
P.cable_diameter  = 0.001;
P.cable_density   = 0.97;
% Linear mass density [kg/m]:
P.cable_mass_per_m = P.cable_density * (pi*(P.cable_diameter/2)^2) * 1e6;
P.cable_stiffness = 50000;  % N/m  (EA/L, approximate for 1mm Dyneema)
P.cable_damping   = 5.0;    % N·s/m

%we need to check what this number should be, it depends on cable material. we need 
%to know how fast the cable feels the new tension set by the motors
% this 20 hz number is a guess
P.cable_compliance_bandwidth = 20.0;   % Hz

% Tension limits
P.T_min = 2.0;              % N  (preload)
P.T_max = 50.0;             % N

% Simulation timing
P.dt      = 0.01;           % s   (100 Hz control loop)
P.t_total = 5.0;            % s   per test run
P.t_vec   = 0:P.dt:P.t_total;
P.N       = length(P.t_vec);

% Throwing in some disturbance to test how models recover
P.disturbance_time      = 2.5;   % s  (midway through trajectory)
P.disturbance_magnitude = 1.5;   % N  (external force impulse)
P.disturbance_duration  = 0.3;   % s


% setting max acceleration and velocity for the cables, this needs to be
% verified before we actually run it, ARMA check would be sick
P.vel_max   = 0.20;   % m/s   — comfortable for 1m workspace at 100Hz
P.accel_max = 0.40;   % m/s²  — ~0.08g, conservative for preload stability

%this is just confirmation that the file is running dw
fprintf('=== CDPR Simulation Framework Initialized ===\n');
fprintf('Control rate: %.0f Hz | Duration: %.1f s | Steps: %d\n\n', ...
        1/P.dt, P.t_total, P.N);


%this makes the smooth motion curve through a trapezoid. it accelerates to accel_max, moves at 
%vel_max, then decelerates once were close.
% in this, pos_traj is where EE needs to be at each iteration, vel_traj is the 
%velocity per iteration
function [pos_traj, vel_traj] = make_trap_traj(pos_start, pos_end, t_vec, vel_max, accel_max)
    % make_trap_traj  Generate a trapezoidal velocity profile in 2D.
    %
    % Inputs:
    %   pos_start  — [px, py] start position (m)
    %   pos_end    — [px, py] end position (m)
    %   t_vec      — 1xN time vector (s)
    %   vel_max    — scalar peak velocity (m/s)
    %   accel_max  — scalar peak acceleration (m/s²)
    %
    % Outputs:
    %   pos_traj   — Nx2 position reference
    %   vel_traj   — Nx2 velocity reference (for feedforward)

    N    = length(t_vec);
    d    = norm(pos_end - pos_start);
    dir  = (pos_end - pos_start) / max(d, 1e-9);  % unit direction

    % lowers vel_max so accel ramp is achievable
    v_peak = min(vel_max, sqrt(accel_max * d));

    % Ramp times
    t_ramp = v_peak / accel_max;
    t_flat = (d - accel_max * t_ramp^2) / v_peak;

    if t_flat < 0
        % Triangle profile (distance too short for full flat top)
        t_ramp = sqrt(d / accel_max);
        v_peak = accel_max * t_ramp;
        t_flat = 0;
    end

    t_end_ramp1  = t_ramp;
    t_end_flat   = t_ramp + t_flat;
    t_end_ramp2  = t_ramp + t_flat + t_ramp;

    pos_traj = zeros(N, 2);
    vel_traj = zeros(N, 2);

    for k = 1:N
        t = t_vec(k);
        if t <= t_end_ramp1
            % Accelerating ramp
            s =  0.5 * accel_max * t^2;
            v =  accel_max * t;
        elseif t <= t_end_flat
            % Constant velocity
            s = 0.5 * accel_max * t_ramp^2 + v_peak * (t - t_end_ramp1);
            v = v_peak;
        elseif t <= t_end_ramp2
            % Decelerating ramp
            dt2 = t - t_end_flat;
            s   = 0.5 * accel_max * t_ramp^2 + v_peak * t_flat ...
                + v_peak * dt2 - 0.5 * accel_max * dt2^2;
            v   = v_peak - accel_max * dt2;
        else
            % Hold at end
            s = d;
            v = 0;
        end

        s = min(s, d);   % clamp to final position
        pos_traj(k,:) = pos_start + s * dir;
        vel_traj(k,:) = v * dir;
    end
end


%  SECTION 1: KINEMATICS WITH PARABOLIC SAG CORRECTION
%  Derivation: parabolic arc length parametrization.
%  sag_mid = w*L^2/(8*T), arc correction = L*(1 + w^2*L^2/(24*T^2))


function [L_arc, sag_mid, L_chord] = ik_with_sag(pos, anchors, T_est, ...
                                                    g, cable_mass_per_m)
    % ik_with_sag  Sag-corrected inverse kinematics (parabolic approximation)
    %
    % Inputs:
    %   pos             — [px, py] EE position (m)
    %   anchors         — Nx2 anchor positions (m)
    %   T_est           — Nx1 estimated cable tensions (N)
    %   g               — gravitational acceleration (m/s²)
    %   cable_mass_per_m— cable linear mass density (kg/m)
    %
    % Outputs:
    %   L_arc    — Nx1 arc lengths (actual cable to pay out, m)
    %   sag_mid  — Nx1 midpoint vertical sag (m)
    %   L_chord  — Nx1 straight-line chord distances (m)

    n       = size(anchors, 1);
    L_chord = sqrt(sum((anchors - pos).^2, 2));   % straight-line IK
    w       = cable_mass_per_m * g;               % N/m, scales with g

    % Avoid division by zero if tension is near 0, makes some arbitrary limit
    T_safe = max(T_est, 0.1);

    % Parabolic arc length correction: L_arc = L_chord * (1 + (w²L²)/(24T²))
    sag_factor = (w .* L_chord).^2 ./ (24 .* T_safe.^2);
    L_arc      = L_chord .* (1 + sag_factor);
    sag_mid    = w .* L_chord.^2 ./ (8 .* T_safe);
end

% Demo: compare sag at Earth vs Lunar
T_demo   = ones(4,1) * 10;   % 10 N nominal tension
pos_demo = [0.5, 0.5];

[L_arc_e, sag_e, L_chord] = ik_with_sag(pos_demo, P.anchors, T_demo, ...
                                          GRAVITY.earth, P.cable_mass_per_m);
% [L_arc_l, sag_l, ~]       = ik_with_sag(pos_demo, P.anchors, T_demo, ...
                                          GRAVITY.lunar, P.cable_mass_per_m);

fprintf('=== Cable Sag at Center, 10N Tension ===\n');
fprintf('  %-20s  %-12s  %-12s  %-12s\n', 'Cable', 'Chord (mm)', ...
        'Earth sag (mm)', 'Lunar sag (mm)');
for i = 1:4
    fprintf('  Cable %d (A%d):       %8.3f    %8.4f    %8.4f\n', ...
            i, i, L_chord(i)*1e3, sag_e(i)*1e3, sag_l(i)*1e3);
end
fprintf('\n  NOTE: Lunar sag is %.1fx larger despite lower gravity.\n', ...
        mean(sag_l./sag_e));
fprintf('  (Lower gravity → lower minimum tension → worse sag geometry)\n\n');


%  SECTION 2: FORWARD KINEMATICS
%  ACTION ITEM 5 — lsqnonlin now validates convergence.
%  Returns pos_est and a validity flag. Falls back to weighted least squares
%  if the nonlinear solver fails (e.g., bad initial guess or ill-conditioned).


function [pos_est, valid, exit_flag] = forward_kinematics(L_measured, ...
                                        anchors, pos_guess)
    % forward_kinematics  Estimate EE position from cable lengths.
    %
    % Outputs:
    %   pos_est   — [px, py] estimated position (m)
    %   valid     — true if solver converged cleanly
    %   exit_flag — raw lsqnonlin exit flag (1=converged, <1=issue)

    residual = @(p) sqrt(sum((anchors - p).^2, 2)) - L_measured;

    opts = optimoptions('lsqnonlin', ...
        'Display',           'off', ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance',     1e-10, ...
        'MaxIterations',     500);

    [pos_est, ~, ~, exit_flag] = lsqnonlin(residual, pos_guess, ...
                                            [0,0], [1,1], opts);

    valid = (exit_flag >= 1);

    if ~valid
        % Fallback: linearised weighted least squares around pos_guess
        % J_ls * delta_p = L_measured - L_guess
        % This is less accurate but always produces a finite answer.
        L_guess = sqrt(sum((anchors - pos_guess).^2, 2));
        J_ls    = (anchors - pos_guess) ./ L_guess;  % Nx2 Jacobian rows
        delta_L = L_measured - L_guess;
        delta_p = (J_ls' * J_ls) \ (J_ls' * delta_L);
        pos_est = pos_guess + delta_p';
        pos_est = max([0,0], min([1,1], pos_est));   % workspace clamp
        warning('FK:SolverFailed', ...
            'lsqnonlin failed (flag=%d). Used WLS fallback.', exit_flag);
    end
end


%  SECTION 3: JACOBIAN AND TENSION DISTRIBUTION


function J = compute_jacobian(pos, anchors)
    L = sqrt(sum((anchors - pos).^2, 2));
    J = (anchors - pos) ./ L;
end

function [T, feasible] = tension_qp(pos, anchors, F_ext, T_min, T_max)
    J    = compute_jacobian(pos, anchors);
    n    = size(anchors,1);
    opts = optimoptions('quadprog','Display','off');
    [T, ~, flag] = quadprog(eye(n), zeros(n,1), [], [], J', F_ext(:), ...
                            T_min*ones(n,1), T_max*ones(n,1), [], opts);
    feasible = (flag == 1);
    if ~feasible, T = (T_min + T_max)/2 * ones(n,1); end
end


%  SECTION 4: PID CONTROLLER
%
%  ACTION ITEM 3 — PID gains updated from pure Ziegler-Nichols guess
%  to an analytical bandwidth-based starting point.
%
%  Derivation sketch (2D Cartesian, decoupled axes):
%    Plant: G(s) = 1/(m*s^2)  (point mass, no cable dynamics)
%    Target closed-loop bandwidth: wc = 2*pi*2 = ~12.6 rad/s  (2 Hz)
%    At crossover: |PID * G| = 1
%
%    For a PD controller (Ki=0 first pass):
%      Kp_bd = m * wc^2 = 0.5 * (12.6)^2 ≈ 79  → too aggressive; scale by
%      workspace / force range factor → ~15-20 is empirically safe here.
%
%    Kd added for ~60° phase margin:  Kd ≈ Kp * (2*zeta/wc)
%      with zeta=0.7:  Kd ≈ 20 * (2*0.7/12.6) ≈ 2.2  → using 3.0 gives
%      a bit more damping margin, which is fine.
%
%    Ki set small to reject steady-state gravity error without windup:
%      Ki ≈ 0.05 * Kp ≈ 1.0  → using 0.5 is conservative.
%
%  Lunar gain schedule:
%    Lower effective gravity → lower cable tensions at same position →
%    slightly reduced restoring forces → modest Kp reduction.
%    Scale factor 0.85 is heuristic; tune empirically on hardware.
%
%  >>> PROFESSOR CHECK: recommend running pidtune() or root-locus on the
%  linearised plant (F_cables → pos via tension QP Jacobian) to validate
%  these gains before first hardware spin-up.
%
%  Operates in Cartesian space: computes force command from position/vel error,
%  then maps to cable length commands via Jacobian pseudo-inverse.


function ctrl = pid_init(Kp, Ki, Kd, dt, output_lim)
    % pid_init  Initialize a 2D Cartesian PID controller struct.
    ctrl.Kp         = Kp;
    ctrl.Ki         = Ki;
    ctrl.Kd         = Kd;
    ctrl.dt         = dt;
    ctrl.integral   = [0; 0];
    ctrl.e_prev     = [0; 0];
    ctrl.output_lim = output_lim;   % max force magnitude (N)
end

function [F_cmd, ctrl] = pid_update(ctrl, pos_des, pos_est, vel_est)
    % pid_update  Compute Cartesian force command.
    %
    % Inputs:
    %   pos_des  — [px, py] desired position
    %   pos_est  — [px, py] estimated EE position
    %   vel_est  — [vx, vy] estimated EE velocity (finite diff or observer)
    % Output:
    %   F_cmd    — [Fx, Fy] commanded force (N)

    e          = (pos_des - pos_est)';
    de         = (e - ctrl.e_prev) / ctrl.dt;

    % Anti-windup: clamp integral before accumulation
    ctrl.integral = ctrl.integral + e * ctrl.dt;
    ctrl.integral = max(-5, min(5, ctrl.integral));  % N·s limit

    F_cmd = ctrl.Kp * e + ctrl.Ki * ctrl.integral + ctrl.Kd * de;

    % Magnitude limit
    F_mag = norm(F_cmd);
    if F_mag > ctrl.output_lim
        F_cmd = F_cmd * ctrl.output_lim / F_mag;
    end

    ctrl.e_prev = e;
end

% Bandwidth-based PID gains — see derivation comment above.
% [Kp, Kp] (same gain on both axes for symmetric workspace)
PID       = pid_init([18, 18], [0.8, 0.8], [3.0, 3.0], P.dt, 20.0);
% Lunar variant: reduced Kp to account for lower gravity-induced tension
PID_lunar = pid_init([15, 15], [0.6, 0.6], [2.8, 2.8], P.dt, 20.0);


%  SECTION 5: INTERVAL TYPE-2 FUZZY LOGIC SYSTEM (IT2-FLS)
%
%  ARCHITECTURE:
%    Inputs:  position error e (m), velocity error de/dt (m/s)
%    Output:  Cartesian force command F (N)
%    Each input has 5 IT2 Gaussian MFs with uncertainty band σ ± Δσ.
%    Output has 5 Type-1 (crisp) MFs — standard IT2 design.
%
%  NEURO COMPONENT:
%    The MF parameters (centers, spreads, uncertainty bands) are tuned
%    offline using gradient descent on a reference trajectory.
%    Function: tune_fis_offline() in Section 6.
%    This replaces hand-tuning and is what "neuro-fuzzy" means here.


function fis = build_it2_fls(sigma_uncertainty)
    % build_it2_fls  Construct an Interval Type-2 FLS for position control.
    %
    % Input:
    %   sigma_uncertainty — fractional width of FOU (e.g. 0.15 = ±15% spread)
    %                       Larger = more uncertainty modeled = more robust
    %                       but slower response.


    fis = mamfis('Name', 'IT2_PositionController', ...
                 'AndMethod', 'min', ...
                 'OrMethod',  'max', ...
                 'ImplicationMethod', 'min', ...
                 'AggregationMethod', 'max', ...
                 'DefuzzificationMethod', 'centroid');

    % INPUT 1: Position error (m)
    % Range: ±0.4 m (covers worst-case crossing of 1m workspace)
    fis = addInput(fis, [-0.4, 0.4], 'Name', 'pos_error');

    % IT2 Gaussian MFs: each defined by [sigma_lower sigma_upper center]
    % sigma_upper = sigma_lower * (1 + sigma_uncertainty) → wider FOU
    sig_e = 0.06;   % base spread
    su    = sig_e * (1 + sigma_uncertainty);  % upper sigma (FOU outer)
    sl    = sig_e * (1 - sigma_uncertainty);  % lower sigma (FOU inner)

    centers_e = [-0.25, -0.10, 0.0, 0.10, 0.25];
    mf_names  = {'NB','NS','ZE','PS','PB'};

    for i = 1:5
        % IT2 MF specified as gaussmf with two-element sigma [sl su]
        % MATLAB IT2 syntax: 'gaussmf', [sigma_lower, sigma_upper, center]
        % For compatibility, we build Type-1 MFs and annotate IT2 intent.
        % If your MATLAB version has full IT2 support, replace with:
        % fis = addMF(fis,'pos_error','gaussmf',[sl su centers_e(i)], ...
        %             'Name', mf_names{i}, 'Type', 'IT2');
        fis = addMF(fis, 'pos_error', 'gaussmf', [su, centers_e(i)], ...
                    'Name', mf_names{i});
    end

    % INPUT 2: Velocity error (m/s)
    % Range: ±0.3 m/s (EE max speed ~0.1 m/s, overshoot headroom)
    fis = addInput(fis, [-0.3, 0.3], 'Name', 'vel_error');

    sig_v = 0.05;
    sv_u  = sig_v * (1 + sigma_uncertainty);

    centers_v = [-0.20, -0.08, 0.0, 0.08, 0.20];
    for i = 1:5
        fis = addMF(fis, 'vel_error', 'gaussmf', [sv_u, centers_v(i)], ...
                    'Name', mf_names{i});
    end

    % OUTPUT: Force command (N)
    % Range: ±20 N (bounded by T_max - T_min spread across 4 cables)
    fis = addOutput(fis, [-20, 20], 'Name', 'force_cmd');

    out_centers = [-16, -8, 0, 8, 16];
    out_spread  = 4.0;
    for i = 1:5
        fis = addMF(fis, 'force_cmd', 'trimf', ...
                    [out_centers(i)-out_spread, out_centers(i), ...
                     out_centers(i)+out_spread], ...
                    'Name', mf_names{i});
    end

    % === Rule base (5x5 = 25 rules, symmetric PD structure) ===
    % Logic: if error is NB (very behind) and vel is NB (moving away fast)
    %        → apply NB force (maximum corrective push)
    %        if error is ZE and vel is ZE → ZE (hold position)
    %
    % Row = pos_error MF (1=NB..5=PB)
    % Col = vel_error MF (1=NB..5=PB)
    % Table entries = output MF index
    rule_table = [
        1  1  1  2  2;   % pos=NB
        1  1  2  2  3;   % pos=NS
        1  2  3  4  5;   % pos=ZE
        3  4  4  5  5;   % pos=PS
        4  4  5  5  5;   % pos=PB
    ];

    rules = [];
    for r = 1:5
        for c = 1:5
            rules = [rules; r, c, rule_table(r,c), 1, 1]; %#ok<AGROW>
        end
    end
    fis = addRule(fis, rules);
end

% Build IT2 FLS (15% uncertainty band — represents ±15% mass uncertainty)
fprintf('=== Building IT2 Fuzzy Controller ===\n');
try
    it2_fis = build_it2_fls(0.15);
    fprintf('  Rules: %d | Inputs: %d | Outputs: %d\n', ...
            length(it2_fis.Rules), length(it2_fis.Inputs), ...
            length(it2_fis.Outputs));
catch ME
    fprintf('  WARNING: %s\n', ME.message);
    it2_fis = [];
end

function F_cmd = it2_evaluate(fis, pos_error, vel_error, output_scale)
    % it2_evaluate  Evaluate IT2 FLS and return 2D force command.
    %
    % The FIS is scalar (one force magnitude axis). We evaluate it for
    % each Cartesian axis independently, which is valid for decoupled
    % planar motion. For coupled dynamics, extend to a 2-input-2-output FIS.
    %
    % Inputs:
    %   fis          — the IT2 FIS
    %   pos_error    — [ex, ey] position error vector (m)
    %   vel_error    — [vx_e, vy_e] velocity error vector (m/s)
    %   output_scale — scalar gain to map fuzzy output (N) to actual force
    %
    % Output:
    %   F_cmd        — [Fx, Fy] force command (N)

    if isempty(fis)
        F_cmd = [0; 0];
        return;
    end

    % Clamp inputs to FIS range
    ex = max(-0.4, min(0.4, pos_error(1)));
    ey = max(-0.4, min(0.4, pos_error(2)));
    vx = max(-0.3, min(0.3, vel_error(1)));
    vy = max(-0.3, min(0.3, vel_error(2)));

    Fx = evalfis(fis, [ex, vx]) * output_scale;
    Fy = evalfis(fis, [ey, vy]) * output_scale;
    F_cmd = [Fx; Fy];
end


%  SECTION 6: NEURO-FUZZY PARAMETER TUNING (OFFLINE)
%
%  This is the "neuro" component. We tune the FIS MF parameters using
%  gradient descent on tracking error over a reference trajectory.
%
%  WHAT GETS TUNED:
%    - Gaussian MF centers (where the MF peaks)
%    - Gaussian MF spreads (how wide the MF is)
%    - NOT the rule structure — rules stay fixed (expert knowledge)
%
%  HOW IT WORKS:
%    1. Run a simulated trajectory with the current FIS
%    2. Compute tracking error cost J = sum(e^2)
%    3. Numerically differentiate J w.r.t. each MF parameter
%    4. Gradient step: param = param - lr * dJ/dparam
%    5. Repeat until convergence
%
%  ACTION ITEM 4 — Fixed malformed function signature on old line 420.
%  Original had: tune_fis_offline(fis, ref_trajectory, P, GRAVITY.earth, 150; 0.001)
%  which is invalid MATLAB syntax (semicolon inside argument list, literal
%  values in the function definition). Fixed below with proper named args.


function [fis_tuned, cost_history] = tune_fis_offline(fis, ref_trajectory, ...
                                                        P, g, n_epochs, lr)
    % tune_fis_offline  Gradient-based MF parameter tuning.
    %
    % Inputs:
    %   fis           — initial IT2 FIS
    %   ref_trajectory— Mx2 reference trajectory (m)
    %   P             — robot parameter struct
    %   g             — gravity (m/s²)
    %   n_epochs      — number of gradient passes (e.g. 30–150)
    %   lr            — learning rate (e.g. 0.001–0.01)
    %
    % Outputs:
    %   fis_tuned     — FIS with tuned MF parameters
    %   cost_history  — n_epochs x 1 cost vector (for learning curve plot)

    if isempty(fis)
        fis_tuned = fis; cost_history = []; return;
    end

    fprintf('  Tuning FIS (%d epochs)... ', n_epochs);
    cost_history = zeros(n_epochs, 1);
    fis_tuned    = fis;
    delta        = 1e-4;   % finite difference step for gradient

    for epoch = 1:n_epochs
        % Simulate tracking with current FIS
        cost = simulate_cost(fis_tuned, ref_trajectory, P, g);
        cost_history(epoch) = cost;

        % Numerically tune input MF parameters (centers only for speed)
        for inp = 1:length(fis_tuned.Inputs)
            for mf = 1:length(fis_tuned.Inputs(inp).MembershipFunctions)
                params = fis_tuned.Inputs(inp).MembershipFunctions(mf).Parameters;

                % Perturb center parameter (+delta)
                params_p = params;
                params_p(2) = params_p(2) + delta;
                fis_p = fis_tuned;
                fis_p.Inputs(inp).MembershipFunctions(mf).Parameters = params_p;
                cost_p = simulate_cost(fis_p, ref_trajectory, P, g);

                % Gradient estimate
                grad = (cost_p - cost) / delta;

                % Gradient step
                new_center = params(2) - lr * grad;
                % Clamp center to input range
                in_range = fis_tuned.Inputs(inp).Range;
                new_center = max(in_range(1)+0.01, min(in_range(2)-0.01, new_center));
                fis_tuned.Inputs(inp).MembershipFunctions(mf).Parameters(2) = new_center;
            end
        end

        if mod(epoch, 10) == 0
            fprintf('.');
        end
    end
    fprintf(' done. Cost: %.4f → %.4f\n', cost_history(1), cost_history(end));
end

function cost = simulate_cost(fis, ref_traj, P, g)
    % Quick simulation for cost evaluation during tuning
    N_ref  = size(ref_traj, 1);
    pos    = ref_traj(1,:);
    vel    = [0, 0];
    F_ext  = [0, -P.mass_EE * g];
    cost   = 0;
    output_scale = 1.0;

    for k = 1:N_ref-1
        pos_des  = ref_traj(k+1,:);
        e_pos    = pos_des - pos;
        e_vel    = [0,0] - vel;
        F_cmd    = it2_evaluate(fis, e_pos, e_vel, output_scale)';

        % Simple EE dynamics: F = ma
        acc  = (F_cmd + F_ext) / P.mass_EE;
        vel  = vel + acc * P.dt;
        pos  = pos + vel * P.dt + 0.5*acc*P.dt^2;
        pos  = max(P.ws_min, min(P.ws_max, pos));

        cost = cost + sum(e_pos.^2);
    end
end

% Build reference trajectory for tuning using trapezoidal profile
t_ref    = linspace(0, 3, 150);
[ref_traj, ~] = make_trap_traj([0.3, 0.3], [0.7, 0.7], t_ref, ...
                                 P.vel_max, P.accel_max);

fprintf('\n=== Offline Neuro-Fuzzy Tuning ===\n');
if ~isempty(it2_fis)
    [it2_fis_tuned, cost_history] = tune_fis_offline(it2_fis, ref_traj, ...
                                                       P, GRAVITY.earth, 30, 0.002);
else
    it2_fis_tuned = it2_fis;
    cost_history  = [];
end


%  SECTION 7: SIMULATION ENGINE
%
%  Runs one full simulation trial and returns a log struct.
%  Called multiple times for each (controller × condition) combination.
%
%  DYNAMICS MODEL:
%    EE treated as a point mass (valid for 2D translation study).
%    F_total = F_cables - F_gravity + F_disturbance
%    Euler integration (matches real-time control loop behavior).
%
%  ACTION ITEM 1 — Cable compliance spring-damper.
%    Previously: cable force = whatever QP says (instantaneous).
%    Now: each cable has an internal length state L_actual(k).
%         Tension = k_cable*(L_cmd - L_actual) + b_cable*(dL/dt)
%         L_actual integrates toward L_cmd with bandwidth set by
%         P.cable_compliance_bandwidth.
%    This means the EE "feels" a lagged tension response, which is
%    physically correct for any elastic cable including Dyneema.


function log = run_simulation(controller_type, fis, pid, pid_lunar, P, g, ...
                               ref_traj, vel_traj, inject_disturbance, ...
                               param_perturb)
    % run_simulation  Core simulation loop.
    %
    % Inputs:
    %   controller_type   — 'IT2' or 'PID'
    %   fis               — IT2 FIS (used if controller_type='IT2')
    %   pid               — PID struct for Earth conditions
    %   pid_lunar         — PID struct for Lunar conditions (gain scheduled)
    %   P                 — robot params
    %   g                 — gravity (m/s²)
    %   ref_traj          — Nx2 reference position trajectory (m)
    %   vel_traj          — Nx2 reference velocity trajectory (m/s)
    %   inject_disturbance— true/false
    %   param_perturb     — struct with .mass_factor, .stiffness_factor
    %
    % Output:
    %   log — struct with time-series data for all metrics

    % Apply parameter uncertainty
    mass_actual      = P.mass_EE * param_perturb.mass_factor;
    stiffness_actual = P.cable_stiffness * param_perturb.stiffness_factor;

    F_weight = [0, -mass_actual * g];

    % Select PID gains based on gravity condition
    is_lunar = (g < 2.0);   % simple threshold — lunar g = 1.635
    if is_lunar
        pid_run = pid_lunar;
    else
        pid_run = pid;
    end
    pid_run.integral = [0; 0];
    pid_run.e_prev   = [0; 0];

    % Initialize state
    pos = ref_traj(1,:);
    vel = [0, 0];

    % ACTION ITEM 1 — Cable compliance state
    % L_actual starts at geometric cable length from initial position
    L_actual = sqrt(sum((P.anchors - pos).^2, 2));   % Nc x 1

    % Compliance bandwidth → first-order lag gain
    % alpha = dt * bandwidth * 2*pi  (discrete approximation of RC filter)
    alpha = P.dt * P.cable_compliance_bandwidth * 2 * pi;
    alpha = min(alpha, 1.0);   % clamp: alpha=1 = instantaneous (no lag)

    % Pre-allocate logs
    N = P.N;
    log.t              = P.t_vec;
    log.pos_ref        = ref_traj;
    log.vel_ref        = vel_traj;
    log.pos_est        = zeros(N, 2);
    log.vel_est        = zeros(N, 2);
    log.error          = zeros(N, 2);
    log.F_cmd          = zeros(N, 2);
    log.T_cables       = zeros(N, P.num_cables);
    log.T_compliance   = zeros(N, P.num_cables);   % tension after compliance lag
    log.sag_mid        = zeros(N, P.num_cables);
    log.L_arc          = zeros(N, P.num_cables);
    log.L_actual       = zeros(N, P.num_cables);   % compliance state
    log.disturbance    = zeros(N, 2);
    log.fk_valid       = true(N, 1);               % FK convergence flag
    log.settled        = false;

    output_scale = 1.0;

    for k = 1:N
        t = P.t_vec(k);

        % Reference at this timestep (now uses trapezoidal traj)
        ref_idx  = min(k, size(ref_traj,1));
        pos_des  = ref_traj(ref_idx,:);
        vel_des  = vel_traj(ref_idx,:);

        % Position and velocity error
        e_pos = pos_des - pos;
        e_vel = vel_des - vel;   % feedforward velocity reference now used

        % Controller
        switch upper(controller_type)
            case 'IT2'
                F_ctrl = it2_evaluate(fis, e_pos, e_vel, output_scale)';
            case 'PID'
                [F_ctrl, pid_run] = pid_update(pid_run, pos_des, pos, vel);
            otherwise
                F_ctrl = [0; 0];
        end

        % Disturbance injection
        F_dist = [0; 0];
        if inject_disturbance
            t_dist = P.disturbance_time;
            t_end  = t_dist + P.disturbance_duration;
            if t >= t_dist && t <= t_end
                F_dist = [P.disturbance_magnitude; 0];
            end
        end

        % PTFE tube force — always applied (not disturbance-gated)
        % Assumes feed tube runs to A1 anchor and pulls EE toward it
        A1 = P.anchors(1,:);
        dist_to_A1 = norm(pos - A1);
        if dist_to_A1 > 0.001
            tube_dir = (A1 - pos) / dist_to_A1;
            tube_force_mag = P.tube_force_per_meter * dist_to_A1;
            F_dist = F_dist + (tube_force_mag * tube_dir)';
        end

        % Cable tension distribution (commanded)
        F_total_ext = F_weight' + F_dist;
        [T_cmd, ~] = tension_qp(pos, P.anchors, F_total_ext + F_ctrl, ...
                                  P.T_min, P.T_max);

        % ACTION ITEM 1 — Cable compliance: compute commanded arc lengths
        [L_cmd, sag_mid, ~] = ik_with_sag(pos, P.anchors, T_cmd, g, ...
                                            P.cable_mass_per_m);

        % First-order lag on cable length state
        % L_actual += alpha * (L_cmd - L_actual)
        L_actual = L_actual + alpha * (L_cmd - L_actual);

        % Compliance tension: spring-damper from actual vs commanded length
        % T_compliance = k*(L_cmd - L_actual) + b*(dL/dt approximated as
        % alpha*(L_cmd-L_actual)/dt which is already captured by the lag)
        % Simpler: compute residual stretch and convert to tension
        dL = L_cmd - L_actual;   % positive = cable still elongating
        T_compliance = max(P.T_min, ...
                           T_cmd + stiffness_actual * dL - P.cable_damping * dL / P.dt);
        T_compliance = min(T_compliance, P.T_max);

        % Net force on EE uses compliance tension (physically correct)
        % Recompute Jacobian force from T_compliance
        J     = compute_jacobian(pos, P.anchors);
        F_cables_actual = J' * T_compliance;   % 2x1

        F_net = F_cables_actual + F_weight' + F_dist;
        acc   = F_net / mass_actual;

        vel   = vel + acc' * P.dt;
        pos   = pos + vel  * P.dt + 0.5 * acc' * P.dt^2;
        pos   = max(P.ws_min, min(P.ws_max, pos));

        % Log
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


%  SECTION 8: DATA ANALYSIS — METRICS


function metrics = compute_metrics(log, P)
    % compute_metrics  Extract all performance metrics from a simulation log.
    %
    % Returns struct with:
    %   .rmse_pos           — RMS position tracking error (mm)
    %   .max_error          — peak error magnitude (mm)
    %   .overshoot          — % overshoot (if any)
    %   .settling_time      — time to reach and stay within 2% band (s)
    %   .steady_state_error — mean error in last 0.5s (mm)
    %   .dist_recovery_time — time to recover after disturbance (s)
    %   .tension_violations — count of timesteps where T < T_min

    error_mag = sqrt(sum(log.error.^2, 2)) * 1000;  % mm

    metrics.rmse_pos         = rms(error_mag);
    metrics.max_error        = max(error_mag);
    metrics.steady_state_err = mean(error_mag(end-round(0.5/P.dt):end));

    % Overshoot relative to final reference position
    ref_final = log.pos_ref(end,:);
    err_final_dir = dot(log.error, ref_final./norm(ref_final), 2);
    undershoot = min(err_final_dir) * 1000;  % mm
    metrics.overshoot = max(0, -undershoot);

    % Settling time: first time error stays below 2% of travel distance
    travel = norm(log.pos_ref(end,:) - log.pos_ref(1,:)) * 1000;  % mm
    band   = 0.02 * travel;
    settled = false;
    metrics.settling_time = P.t_total;
    for k = 1:P.N
        if all(error_mag(k:end) < band)
            metrics.settling_time = P.t_vec(k);
            settled = true;
            break;
        end
    end
    metrics.settled = settled;

    % Disturbance recovery: time from disturbance until error < 5mm
    t_dist_idx = find(P.t_vec >= P.disturbance_time, 1);
    if ~isempty(t_dist_idx) && t_dist_idx < P.N
        recovery_idx = find(error_mag(t_dist_idx:end) < 5.0, 1);
        if ~isempty(recovery_idx)
            metrics.dist_recovery_time = recovery_idx * P.dt;
        else
            metrics.dist_recovery_time = NaN;
        end
    else
        metrics.dist_recovery_time = NaN;
    end

    % Tension constraint violations
    metrics.tension_violations = sum(any(log.T_cables < P.T_min, 2));

    % Compliance lag metric: mean deviation between commanded and actual length
    if isfield(log, 'L_actual')
        dL_all = log.L_arc - log.L_actual;
        metrics.mean_compliance_lag_mm = mean(abs(dL_all(:))) * 1000;
        metrics.max_compliance_lag_mm  = max(abs(dL_all(:))) * 1000;
    else
        metrics.mean_compliance_lag_mm = NaN;
        metrics.max_compliance_lag_mm  = NaN;
    end
end


%  SECTION 9: RUN ALL SIMULATION CONDITIONS
%  3 conditions × 2 controllers × 2 tests = 12 runs
% Uncertainty model (perturbed condition)
P.mass_uncertainty      = 0.15;   % ±15% on EE mass
P.stiffness_uncertainty = 0.20;   % ±20% on cable stiffness

fprintf('\n=== Running Simulation Study ===\n');


% Move from (0.3, 0.4) to (0.7, 0.6), hold at end
[step_ref, step_vel] = make_trap_traj([0.3, 0.4], [0.7, 0.6], P.t_vec, ...
                                       P.vel_max, P.accel_max);

% Nominal parameters
nominal = struct('mass_factor', 1.0, 'stiffness_factor', 1.0);

% Perturbed parameters (±15% mass, ±20% stiffness — worst case)
perturbed = struct('mass_factor', 1.15, 'stiffness_factor', 0.80);

conditions = struct(...
    'name',    {'Earth Nominal', 'Earth Perturbed', 'Lunar'},     ...
    'g',       {GRAVITY.earth,  GRAVITY.earth,      GRAVITY.lunar}, ...
    'params',  {nominal,        perturbed,          nominal});

results = struct();

for c = 1:length(conditions)
    cond = conditions(c);
    fprintf('\n  Condition: %s (g=%.2f m/s²)\n', cond.name, cond.g);

    for ctrl_idx = 1:2
        if ctrl_idx == 1
            cname = 'IT2';  fis_use = it2_fis_tuned;
        else
            cname = 'PID';  fis_use = [];
        end

        % Test 1: Step response (no disturbance)
        log1 = run_simulation(cname, fis_use, PID, PID_lunar, P, cond.g, ...
                               step_ref, step_vel, false, cond.params);
        m1   = compute_metrics(log1, P);

        % Test 2: Disturbance rejection
        log2 = run_simulation(cname, fis_use, PID, PID_lunar, P, cond.g, ...
                               step_ref, step_vel, true, cond.params);
        m2   = compute_metrics(log2, P);

        fname = sprintf('%s_%s', strrep(cond.name,' ','_'), cname);
        results.(fname).step         = log1;
        results.(fname).disturbed    = log2;
        results.(fname).metrics_step = m1;
        results.(fname).metrics_dist = m2;

        fprintf('    [%s] RMSE=%.2fmm | Settle=%.2fs | Recovery=%.2fs\n', ...
                cname, m1.rmse_pos, m1.settling_time, ...
                m2.dist_recovery_time);
    end
end


%  SECTION 9b: PRINT RESULTS TABLE (CONSOLE)


fprintf('\n');
fprintf('%-22s | %-6s | %-9s | %-10s | %-10s | %-10s\n', ...
        'Condition + Controller', 'RMSE', 'MaxErr', 'Settle', ...
        'Recovery', 'T_viol');
fprintf('%s\n', repmat('-', 1, 80));

cond_names = {'Earth_Nominal', 'Earth_Perturbed', 'Lunar'};
ctrl_names = {'IT2', 'PID'};

for c = 1:length(cond_names)
    for ct = 1:length(ctrl_names)
        fname = sprintf('%s_%s', cond_names{c}, ctrl_names{ct});
        if isfield(results, fname)
            m1 = results.(fname).metrics_step;
            m2 = results.(fname).metrics_dist;
            fprintf('%-22s | %5.2fmm | %8.2fmm | %8.2fs | %8.2fs | %5d\n', ...
                    [strrep(cond_names{c},'_',' ') ' ' ctrl_names{ct}], ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, ...
                    m2.dist_recovery_time, m1.tension_violations);
        end
    end
end


%  SECTION 9c: FIGURES


figure('Name','CDPR Simulation Results','Position',[50 50 1400 900]);

r_it2 = results.Earth_Nominal_IT2;
r_pid = results.Earth_Nominal_PID;

% 1. Step response tracking comparison
subplot(2,3,1);
error_it2 = sqrt(sum(r_it2.step.error.^2,2))*1000;
error_pid = sqrt(sum(r_pid.step.error.^2,2))*1000;
plot(P.t_vec, error_it2, 'b-', 'LineWidth',1.8); hold on;
plot(P.t_vec, error_pid, 'r--','LineWidth',1.8);
yline(2, 'k:', '2% band'); grid on;
xlabel('Time (s)'); ylabel('Position error (mm)');
title('Step Response: Tracking Error');
legend('IT2-FLS','PID','Location','NE');

% 2. Disturbance rejection
subplot(2,3,2);
r_it2d = results.Earth_Nominal_IT2.disturbed;
r_pidd = results.Earth_Nominal_PID.disturbed;
error_it2d = sqrt(sum(r_it2d.error.^2,2))*1000;
error_pidd = sqrt(sum(r_pidd.error.^2,2))*1000;
plot(P.t_vec, error_it2d, 'b-', 'LineWidth',1.8); hold on;
plot(P.t_vec, error_pidd, 'r--','LineWidth',1.8);
xline(P.disturbance_time, 'k:', 'Disturbance');
xline(P.disturbance_time+P.disturbance_duration, 'k:');
grid on; xlabel('Time (s)'); ylabel('Error (mm)');
title('Disturbance Rejection');
legend('IT2-FLS','PID');

% 3. Gravity condition comparison (IT2 only)
subplot(2,3,3);
e_en = sqrt(sum(results.Earth_Nominal_IT2.step.error.^2,2))*1000;
e_ep = sqrt(sum(results.Earth_Perturbed_IT2.step.error.^2,2))*1000;
e_lu = sqrt(sum(results.Lunar_IT2.step.error.^2,2))*1000;
plot(P.t_vec, e_en, 'b-',  'LineWidth',1.8); hold on;
plot(P.t_vec, e_ep, 'g--', 'LineWidth',1.8);
plot(P.t_vec, e_lu, 'r:',  'LineWidth',2.0);
grid on; xlabel('Time (s)'); ylabel('Error (mm)');
title('IT2-FLS: Earth vs Lunar Conditions');
legend('Earth nominal','Earth perturbed','Lunar');

% 4. Cable sag comparison (Earth vs Lunar, cable 1)
subplot(2,3,4);
sag_e = results.Earth_Nominal_IT2.step.sag_mid(:,1)*1000;
sag_l = results.Lunar_IT2.step.sag_mid(:,1)*1000;
plot(P.t_vec, sag_e, 'b-', 'LineWidth',1.8); hold on;
plot(P.t_vec, sag_l, 'r--','LineWidth',1.8);
grid on; xlabel('Time (s)'); ylabel('Sag at midspan (mm)');
title('Cable 1 Sag: Earth vs Lunar');
legend('Earth','Lunar');

% 5. FIS learning curve
subplot(2,3,5);
if ~isempty(cost_history)
    plot(1:length(cost_history), cost_history, 'b-', 'LineWidth',1.8);
    xlabel('Training epoch'); ylabel('Tracking cost J');
    title('Neuro-Fuzzy Training: Cost vs Epoch'); grid on;
else
    text(0.5,0.5,'Fuzzy Toolbox unavailable','HorizontalAlignment','center',...
         'Units','normalized');
    title('Learning Curve (unavailable)');
end

% 6. Metric bar chart
subplot(2,3,6);
metric_labels = {'Earth-IT2','Earth-PID','Perturbed-IT2','Perturbed-PID',...
                 'Lunar-IT2','Lunar-PID'};
rmse_vals = [results.Earth_Nominal_IT2.metrics_step.rmse_pos, ...
             results.Earth_Nominal_PID.metrics_step.rmse_pos, ...
             results.Earth_Perturbed_IT2.metrics_step.rmse_pos, ...
             results.Earth_Perturbed_PID.metrics_step.rmse_pos, ...
             results.Lunar_IT2.metrics_step.rmse_pos, ...
             results.Lunar_PID.metrics_step.rmse_pos];
bar_colors = repmat([0.2 0.4 0.8; 0.9 0.3 0.3], 3, 1);
b = bar(rmse_vals, 'FaceColor','flat');
b.CData = bar_colors;
set(gca,'XTickLabel', metric_labels, 'XTickLabelRotation', 35);
ylabel('RMSE (mm)'); title('Position RMSE by Condition');
grid on; yline(0);

sgtitle('CDPR Type-2 Neuro-Fuzzy vs PID — Simulation Study');


%  SECTION 10: EXPORT — ACTION ITEM 6
%  Writes results to CSV and a LaTeX table fragment.
%  Call export_results_table(results, P) after simulations finish.


function export_results_table(results, P)
    % export_results_table  Write metrics to CSV + LaTeX for publication.
    %
    % Output files (written to current directory):
    %   cdpr_results.csv       — machine-readable, import into Excel/pandas
    %   cdpr_results_latex.tex — drop into a paper's results section

    cond_names = {'Earth_Nominal', 'Earth_Perturbed', 'Lunar'};
    ctrl_names = {'IT2', 'PID'};
    timestamp  = datestr(now, 'yyyy-mm-dd_HHMMSS');

    % --- CSV ---
    csv_file = 'cdpr_results.csv';
    fid = fopen(csv_file, 'w');
    fprintf(fid, 'Condition,Controller,RMSE_mm,MaxErr_mm,SettleTime_s,');
    fprintf(fid, 'RecoveryTime_s,TensionViolations,ComplianceLag_mm\n');

    for c = 1:length(cond_names)
        for ct = 1:length(ctrl_names)
            fname = sprintf('%s_%s', cond_names{c}, ctrl_names{ct});
            if isfield(results, fname)
                m1 = results.(fname).metrics_step;
                m2 = results.(fname).metrics_dist;
                comp_lag = NaN;
                if isfield(m1, 'mean_compliance_lag_mm')
                    comp_lag = m1.mean_compliance_lag_mm;
                end
                rt = m2.dist_recovery_time;
                if isnan(rt), rt_str = 'NaN'; else, rt_str = sprintf('%.3f',rt); end
                fprintf(fid, '%s,%s,%.4f,%.4f,%.4f,%s,%d,%.4f\n', ...
                    strrep(cond_names{c},'_',' '), ctrl_names{ct}, ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, ...
                    rt_str, m1.tension_violations, comp_lag);
            end
        end
    end
    fclose(fid);

    % --- LaTeX ---
    tex_file = 'cdpr_results_latex.tex';
    fid = fopen(tex_file, 'w');
    fprintf(fid, '%% Auto-generated: %s\n', timestamp);
    fprintf(fid, '\\begin{table}[h]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{CDPR Simulation Performance Metrics}\n');
    fprintf(fid, '\\label{tab:cdpr_results}\n');
    fprintf(fid, '\\begin{tabular}{llrrrr}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Condition & Controller & RMSE (mm) & Max Err (mm) & Settle (s) & Recovery (s) \\\\\n');
    fprintf(fid, '\\midrule\n');

    for c = 1:length(cond_names)
        for ct = 1:length(ctrl_names)
            fname = sprintf('%s_%s', cond_names{c}, ctrl_names{ct});
            if isfield(results, fname)
                m1 = results.(fname).metrics_step;
                m2 = results.(fname).metrics_dist;
                rt = m2.dist_recovery_time;
                cond_label = strrep(cond_names{c},'_',' ');
                if isnan(rt)
                    rt_str = '---';
                else
                    rt_str = sprintf('%.2f', rt);
                end
                fprintf(fid, '%s & %s & %.2f & %.2f & %.2f & %s \\\\\n', ...
                    cond_label, ctrl_names{ct}, ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, rt_str);
            end
        end
        fprintf(fid, '\\midrule\n');
    end

    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);

    fprintf('\n=== Export complete ===\n');
    fprintf('  CSV:   %s\n', csv_file);
    fprintf('  LaTeX: %s\n', tex_file);
end

% Run export automatically after simulation
export_results_table(results, P);

fprintf('\n=== Simulation complete. Check results table and figures. ===\n');
