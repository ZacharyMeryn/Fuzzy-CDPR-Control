
% impactful asf
% Type-2 Neuro-Fuzzy vs PID Simulation Study — EARTH CONDITIONS ONLY
% Cable-Driven Parallel Robot — 1m x 1m
%
% PURPOSE OF THIS FILE:
%   Compare IT2-FLS vs PID on just earths gravity, no cable sag
%   Two conditions: nominal parameters, and perturbed parameters (simulates
%   slight model mismatch, good for testing). Two tests per condition: step
%  response and disturbance recalibration
%   Total: 2 conditions x 2 controllers x 2 tests = 8 runs.
%
% sequel file is on the github
%
% BEFORE RUNNING ANYTHING, FILL IN MEASURED VALUES:
%   P.mass_EE          — weigh the end effector (kg)
%   P.cable_diameter   — measure with calipers (m)
%   P.cable_mass_per_m — from Dyneema spec sheet (kg/m)
%   P.cable_stiffness  — axial stiffness EA/L from spec sheet (N/m)
%   P.T_min / P.T_max  — from motor rated torque and cable break strength
%
% BEFORE RUNNING, DO THE FOLLOWING:
%   1. Reorder function definitions to the bottom of the file
%   2. Confirm Optimization Toolbox and Fuzzy Logic Toolbox are installed (run: ver)
%   3. Run smoke test with placeholder values
%
% REQUIRED TOOLBOXES:
%   - Optimization Toolbox  (quadprog, lsqnonlin)
%   - Fuzzy Logic Toolbox   (R2019b+ for Type-2 support)
%
% CABLE SAG MODEL: Parabolic approximation
%   sag_mid = w*L^2/(8*T), arc correction = L*(1 + w^2*L^2/(24*T^2))
%   where w = cable_mass_per_m * g

clear; clc; close all;


%  SECTION 0: PARAMETERS

GRAVITY.earth = 9.81;   % m/s^2

% Workspace geometry — 1m x 1m, motors at corners
P.anchors = [0.0, 1.0;      % A1 top left
             1.0, 1.0;      % A2 top right
             1.0, 0.0;      % A3 bottom right
             0.0, 0.0];     % A4 bottom left
P.num_cables = 4;
P.ws_min = [0.05, 0.05];
P.ws_max = [0.95, 0.95];

% End-effector — placeholder, measure before running for real
P.mass_EE    = 0.5;      % kg
P.inertia_EE = 0.001;   % kg*m^2

% PTFE tube force pulling EE toward A1 anchor
% units N/m — measure when assembled
P.tube_force_per_meter = 0.3;

% Cable properties — placeholder, get from spec sheet
P.cable_diameter  = 0.001;   % m
P.cable_density   = 0.97;    % kg/m^3
P.cable_mass_per_m = P.cable_density * (pi*(P.cable_diameter/2)^2) * 1e6;
P.cable_stiffness = 50000;   % N/m (EA/L, approximate for 1mm Dyneema)
P.cable_damping   = 5.0;     % N*s/m

% How fast cables feel a tension change — needs prof check
% 20 Hz is a guess based on Dyneema stiffness characteristics
P.cable_compliance_bandwidth = 20.0;   % Hz

% Tension limits — derive from motor torque and cable break strength
P.T_min = 2.0;    % N (preload — keeps all cables taut)
P.T_max = 50.0;   % N

% Timing
P.dt      = 0.01;           % s (100 Hz control loop)
P.t_total = 5.0;            % s per run
P.t_vec   = 0:P.dt:P.t_total;
P.N       = length(P.t_vec);

% Disturbance — a sideways shove at 2.5s lasting 0.3s
P.disturbance_time      = 2.5;   % s
P.disturbance_magnitude = 1.5;   % N
P.disturbance_duration  = 0.3;   % s

% Trajectory speed limits — verify before hardware runs
% >>> PROFESSOR CHECK: confirm vel_max and accel_max before any motor input
P.vel_max   = 0.20;   % m/s
P.accel_max = 0.40;   % m/s^2

fprintf('=== CDPR Earth Simulation Initialized ===\n');
fprintf('Control rate: %.0f Hz | Duration: %.1f s | Steps: %d\n\n', ...
        1/P.dt, P.t_total, P.N);


%  Prelim part 2: smooth operator
%  smooth movement, stops teleporting between start and end points
%  Accelerates to vel_max, cruises, decelerates. If distance is short,
%  uses a triangle profile, no cruising altitude

function [pos_traj, vel_traj] = make_trap_traj(pos_start, pos_end, t_vec, vel_max, accel_max)
    N    = length(t_vec);
    d    = norm(pos_end - pos_start);
    dir  = (pos_end - pos_start) / max(d, 1e-9);

    v_peak = min(vel_max, sqrt(accel_max * d));
    t_ramp = v_peak / accel_max;
    t_flat = (d - accel_max * t_ramp^2) / v_peak;

    if t_flat < 0
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
            s = d;
            v = 0;
        end
        s = min(s, d);
        pos_traj(k,:) = pos_start + s * dir;
        vel_traj(k,:) = v * dir;
    end
end


%  SECTION 1: INVERSE KINEMATICS WITH PARABOLIC SAG CORRECTION
%  Given where the EE is, how long does each cable need to be
%  adds straight line distance and the sag correction

function [L_arc, sag_mid, L_chord] = ik_with_sag(pos, anchors, T_est, g, cable_mass_per_m)
    L_chord    = sqrt(sum((anchors - pos).^2, 2));
    w          = cable_mass_per_m * g;
    T_safe     = max(T_est, 0.1);
    sag_factor = (w .* L_chord).^2 ./ (24 .* T_safe.^2);
    L_arc      = L_chord .* (1 + sag_factor);
    sag_mid    = w .* L_chord.^2 ./ (8 .* T_safe);
end

% dummy test to check if everythings working, just puts sample strain on
% all cables to see if they work ok
T_demo   = ones(4,1) * 10;
pos_demo = [0.5, 0.5];
[~, sag_e, L_chord] = ik_with_sag(pos_demo, P.anchors, T_demo, ...
                                    GRAVITY.earth, P.cable_mass_per_m);
fprintf('=== Cable Sag at Center (Earth, 10N) ===\n');
for i = 1:4
    fprintf('  Cable %d: chord=%.1fmm  sag=%.4fmm\n', ...
            i, L_chord(i)*1e3, sag_e(i)*1e3);
end
fprintf('\n');


%  SECTION 2: FORWARD KINEMATICS
%  Given cable lengths, where is the EE
%  Uses lsqnonlin, then weighted least squares if that fails

function [pos_est, valid, exit_flag] = forward_kinematics(L_measured, anchors, pos_guess)
    residual = @(p) sqrt(sum((anchors - p).^2, 2)) - L_measured;
    opts = optimoptions('lsqnonlin', 'Display', 'off', ...
        'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, 'MaxIterations', 500);
    [pos_est, ~, ~, exit_flag] = lsqnonlin(residual, pos_guess, [0,0], [1,1], opts);
    valid = (exit_flag >= 1);
    if ~valid
        L_guess = sqrt(sum((anchors - pos_guess).^2, 2));
        J_ls    = (anchors - pos_guess) ./ L_guess;
        delta_L = L_measured - L_guess;
        delta_p = (J_ls' * J_ls) \ (J_ls' * delta_L);
        pos_est = max([0,0], min([1,1], pos_guess + delta_p'));
        warning('FK:SolverFailed', 'lsqnonlin failed (flag=%d). Used WLS fallback.', exit_flag);
    end
end


%  SECTION 3: JACOBIAN AND TENSION DISTRIBUTION
%  Jacobian maps cable tensions to EE force.
%  QP finds the minimum-norm tension set that produces the commanded force
%  while keeping all tensions inside [T_min, T_max].

function J = compute_jacobian(pos, anchors)
    L = sqrt(sum((anchors - pos).^2, 2));
    J = (anchors - pos) ./ L;
end

function [T, feasible] = tension_qp(pos, anchors, F_ext, T_min, T_max)
    J    = compute_jacobian(pos, anchors);
    n    = size(anchors, 1);
    opts = optimoptions('quadprog', 'Display', 'off');
    [T, ~, flag] = quadprog(eye(n), zeros(n,1), [], [], J', F_ext(:), ...
                            T_min*ones(n,1), T_max*ones(n,1), [], opts);
    feasible = (flag == 1);
    if ~feasible, T = (T_min + T_max)/2 * ones(n,1); end
end


%  SECTION 4: PID CONTROLLER
%  Bandwidth-based gain derivation (2 Hz target bandwidth, 0.5 kg EE):
%    Kp ~ m * wc^2 scaled to force range ~ 18
%    Kd ~ Kp * 2*zeta/wc with zeta=0.7 ~ 3.0
%    Ki small to reject gravity offset without windup ~ 0.8
%
%  ask ARMA to validate gains with pidtune() or root-locus before
%  running, especially with hardware

function ctrl = pid_init(Kp, Ki, Kd, dt, output_lim)
    ctrl.Kp = Kp;  ctrl.Ki = Ki;  ctrl.Kd = Kd;
    ctrl.dt = dt;
    ctrl.integral  = [0; 0];
    ctrl.e_prev    = [0; 0];
    ctrl.output_lim = output_lim;
end

function [F_cmd, ctrl] = pid_update(ctrl, pos_des, pos_est, vel_est)
    e             = (pos_des - pos_est)';
    de            = (e - ctrl.e_prev) / ctrl.dt;
    ctrl.integral = ctrl.integral + e * ctrl.dt;
    ctrl.integral = max(-5, min(5, ctrl.integral));
    F_cmd = ctrl.Kp * e + ctrl.Ki * ctrl.integral + ctrl.Kd * de;
    F_mag = norm(F_cmd);
    if F_mag > ctrl.output_lim
        F_cmd = F_cmd * ctrl.output_lim / F_mag;
    end
    ctrl.e_prev = e;
end

PID = pid_init([18, 18], [0.8, 0.8], [3.0, 3.0], P.dt, 20.0);


%  SECTION 5: IT2 FUZZY LOGIC SYSTEM
%  5 MFs per input (NB NS ZE PS PB), 25-rule table.
%  sigma_uncertainty sets the width of the F.o.U
%  larger = more robust to noise, slower response.
% the 5 inputs are NB (negative big thats far behind)
%                  NS (negative small its less behind)
%                  ZE (ideally on target)
%                  PS (positive small, u get it)
%                  PB (positive big)

function fis = build_it2_fls(sigma_uncertainty)
    fis = mamfis('Name', 'IT2_PositionController', ...
                 'AndMethod', 'min', 'OrMethod', 'max', ...
                 'ImplicationMethod', 'min', 'AggregationMethod', 'max', ...
                 'DefuzzificationMethod', 'centroid');
% adds position error, -0.4,0.4 are the bounds of how far that error
% would be, fair enough since it covers sm of the workspace. 
    fis = addInput(fis, [-0.4, 0.4], 'Name', 'pos_error');
    sig_e = 0.06;     % gaussian base width
    su    = sig_e * (1 + sigma_uncertainty);
    centers_e = [-0.25, -0.10, 0.0, 0.10, 0.25]; % centers for the 5 functions
    mf_names  = {'NB','NS','ZE','PS','PB'};
    for i = 1:5
        fis = addMF(fis, 'pos_error', 'gaussmf', [su, centers_e(i)], 'Name', mf_names{i});
    end %that just added the gaussian bell curves for the MFs

    fis = addInput(fis, [-0.3, 0.3], 'Name', 'vel_error');
    sig_v = 0.05;
    sv_u  = sig_v * (1 + sigma_uncertainty);
    centers_v = [-0.20, -0.08, 0.0, 0.08, 0.20];
    for i = 1:5
        fis = addMF(fis, 'vel_error', 'gaussmf', [sv_u, centers_v(i)], 'Name', mf_names{i});
    end
% this controls the output, its a force interval from -20 to 20 N. its a
% triangle instead of a gaussian bell curve bc mamdani systems prefer
% computationally cheap formats. 
    fis = addOutput(fis, [-20, 20], 'Name', 'force_cmd');
    out_centers = [-16, -8, 0, 8, 16];
    for i = 1:5
        fis = addMF(fis, 'force_cmd', 'trimf', ...
                    [out_centers(i)-4, out_centers(i), out_centers(i)+4], 'Name', mf_names{i});
    end
% this looks weird at first, TLDR more left means far behind position wise,
% more up means far behind velocity wise.
    rule_table = [1 1 1 2 2; 
                  1 1 2 2 3; 
                  1 2 3 4 5; 
                  3 4 4 5 5; 
                  4 4 5 5 5];
    rules = [];
    for r = 1:5
        for c = 1:5
            rules = [rules; r, c, rule_table(r,c), 1, 1]; %#ok<AGROW>
        end
    end
    fis = addRule(fis, rules);
end

fprintf('=== Building IT2 Fuzzy Controller ===\n');
try
    it2_fis = build_it2_fls(0.15);
    fprintf('  Rules: %d | Inputs: %d | Outputs: %d\n', ...
            length(it2_fis.Rules), length(it2_fis.Inputs), length(it2_fis.Outputs));
catch ME
    fprintf('  WARNING: %s\n', ME.message);
    it2_fis = [];
end

function F_cmd = it2_evaluate(fis, pos_error, vel_error, output_scale)
    if isempty(fis), F_cmd = [0; 0]; return; end
    ex = max(-0.4, min(0.4, pos_error(1)));
    ey = max(-0.4, min(0.4, pos_error(2)));
    vx = max(-0.3, min(0.3, vel_error(1)));
    vy = max(-0.3, min(0.3, vel_error(2)));
    Fx = evalfis(fis, [ex, vx]) * output_scale;
    Fy = evalfis(fis, [ey, vy]) * output_scale;
    F_cmd = [Fx; Fy];
end


%  SECTION 6: NEURO-FUZZY OFFLINE TUNING
%  Gradient descent on MF centers to minimize tracking cost.
%  pushes each center, checks cost, and then re adjusts as needed

function [fis_tuned, cost_history] = tune_fis_offline(fis, ref_trajectory, P, g, n_epochs, lr)
    if isempty(fis), fis_tuned = fis; cost_history = []; return; end % checks if you have the right library in
    fprintf('  Tuning FIS (%d epochs)... ', n_epochs);
    cost_history = zeros(n_epochs, 1); %creates array to hold training data
    fis_tuned    = fis; % original mamfis copy
    delta        = 1e-4;
    for epoch = 1:n_epochs % epoch -> pass through all MFs, produces single squared tracking error
        cost = simulate_cost(fis_tuned, ref_trajectory, P, g);
        cost_history(epoch) = cost;
        for inp = 1:length(fis_tuned.Inputs)
            for mf = 1:length(fis_tuned.Inputs(inp).MembershipFunctions)
                params   = fis_tuned.Inputs(inp).MembershipFunctions(mf).Parameters;
                params_p = params;
                params_p(2) = params_p(2) + delta;
                fis_p    = fis_tuned;
                fis_p.Inputs(inp).MembershipFunctions(mf).Parameters = params_p;
                cost_p   = simulate_cost(fis_p, ref_trajectory, P, g);
                grad     = (cost_p - cost) / delta;
                new_center = params(2) - lr * grad;
                in_range   = fis_tuned.Inputs(inp).Range;
                new_center = max(in_range(1)+0.01, min(in_range(2)-0.01, new_center));
                fis_tuned.Inputs(inp).MembershipFunctions(mf).Parameters(2) = new_center;
            end
        end
        if mod(epoch, 10) == 0, fprintf('.'); end
    end
    fprintf(' done. Cost: %.4f -> %.4f\n', cost_history(1), cost_history(end));
end
% everything above this just grabs each MF centroid
function cost = simulate_cost(fis, ref_traj, P, g)
    N_ref = size(ref_traj, 1);
    pos   = ref_traj(1,:);
    vel   = [0, 0];
    F_ext = [0, -P.mass_EE * g];
    cost  = 0;
    for k = 1:N_ref-1
        pos_des = ref_traj(k+1,:);
        e_pos   = pos_des - pos;
        e_vel   = [0,0] - vel;
        F_cmd   = it2_evaluate(fis, e_pos, e_vel, 1.0)';
        acc     = (F_cmd + F_ext) / P.mass_EE;
        vel     = vel + acc * P.dt;
        pos     = pos + vel * P.dt + 0.5*acc*P.dt^2;
        pos     = max(P.ws_min, min(P.ws_max, pos));
        cost    = cost + sum(e_pos.^2);
    end
end

t_ref = linspace(0, 3, 150);
[ref_traj_tune, ~] = make_trap_traj([0.3, 0.3], [0.7, 0.7], t_ref, P.vel_max, P.accel_max);

fprintf('\n=== Offline Neuro-Fuzzy Tuning (Earth) ===\n');
if ~isempty(it2_fis)
    [it2_fis_tuned, cost_history] = tune_fis_offline(it2_fis, ref_traj_tune, ...
                                                       P, GRAVITY.earth, 30, 0.002);
else
    it2_fis_tuned = it2_fis;
    cost_history  = [];
end


%  SECTION 7: SIMULATION ENGINE
%  One call per trial. Runs the full timestep loop, logs everything.
%  Cable compliance: L_actual lags behind L_cmd at the bandwidth set above.
%  The EE feels the lagged tension, not the instantaneous commanded tension.

function log = run_simulation(controller_type, fis, pid, P, g, ...
                               ref_traj, vel_traj, inject_disturbance, param_perturb)
    mass_actual      = P.mass_EE * param_perturb.mass_factor;
    stiffness_actual = P.cable_stiffness * param_perturb.stiffness_factor;
    F_weight         = [0, -mass_actual * g];

    pid_run          = pid;
    pid_run.integral = [0; 0];
    pid_run.e_prev   = [0; 0];

    pos      = ref_traj(1,:);
    vel      = [0, 0];
    L_actual = sqrt(sum((P.anchors - pos).^2, 2));
    alpha    = min(P.dt * P.cable_compliance_bandwidth * 2 * pi, 1.0);

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

        switch upper(controller_type)
            case 'IT2'
                F_ctrl = it2_evaluate(fis, e_pos, e_vel, 1.0)';
            case 'PID'
                [F_ctrl, pid_run] = pid_update(pid_run, pos_des, pos, vel);
            otherwise
                F_ctrl = [0; 0];
        end

        F_dist = [0; 0];
        if inject_disturbance
            if t >= P.disturbance_time && t <= P.disturbance_time + P.disturbance_duration
                F_dist = [P.disturbance_magnitude; 0];
            end
        end

        % PTFE tube force toward A1
        A1 = P.anchors(1,:);
        dist_to_A1 = norm(pos - A1);
        if dist_to_A1 > 0.001
            tube_dir = (A1 - pos) / dist_to_A1;
            F_dist = F_dist + (P.tube_force_per_meter * dist_to_A1 * tube_dir)';
        end

        F_total_ext = F_weight' + F_dist;
        [T_cmd, ~]  = tension_qp(pos, P.anchors, F_total_ext + F_ctrl, P.T_min, P.T_max);

        [L_cmd, sag_mid, ~] = ik_with_sag(pos, P.anchors, T_cmd, g, P.cable_mass_per_m);

        L_actual = L_actual + alpha * (L_cmd - L_actual);
        dL       = L_cmd - L_actual;
        T_compliance = max(P.T_min, min(P.T_max, ...
                           T_cmd + stiffness_actual * dL - P.cable_damping * dL / P.dt));

        J               = compute_jacobian(pos, P.anchors);
        F_cables_actual = J' * T_compliance;
        F_net           = F_cables_actual + F_weight' + F_dist;
        acc             = F_net / mass_actual;

        vel = vel + acc' * P.dt;
        pos = pos + vel * P.dt + 0.5 * acc' * P.dt^2;
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


%  SECTION 8: METRICS

function metrics = compute_metrics(log, P)
    error_mag = sqrt(sum(log.error.^2, 2)) * 1000;   % mm
    metrics.rmse_pos         = rms(error_mag);
    metrics.max_error        = max(error_mag);
    metrics.steady_state_err = mean(error_mag(end-round(0.5/P.dt):end));

    ref_final     = log.pos_ref(end,:);
    err_final_dir = dot(log.error, ref_final./norm(ref_final), 2);
    metrics.overshoot = max(0, -min(err_final_dir)*1000);

    travel = norm(log.pos_ref(end,:) - log.pos_ref(1,:)) * 1000;
    band   = 0.02 * travel;
    metrics.settling_time = P.t_total;
    for k = 1:P.N
        if all(error_mag(k:end) < band)
            metrics.settling_time = P.t_vec(k);
            break;
        end
    end

    t_dist_idx = find(P.t_vec >= P.disturbance_time, 1);
    if ~isempty(t_dist_idx) && t_dist_idx < P.N
        recovery_idx = find(error_mag(t_dist_idx:end) < 5.0, 1);
        metrics.dist_recovery_time = recovery_idx * P.dt;
        if isempty(recovery_idx), metrics.dist_recovery_time = NaN; end
    else
        metrics.dist_recovery_time = NaN;
    end

    metrics.tension_violations = sum(any(log.T_cables < P.T_min, 2));

    if isfield(log, 'L_actual')
        dL_all = log.L_arc - log.L_actual;
        metrics.mean_compliance_lag_mm = mean(abs(dL_all(:))) * 1000;
        metrics.max_compliance_lag_mm  = max(abs(dL_all(:))) * 1000;
    else
        metrics.mean_compliance_lag_mm = NaN;
        metrics.max_compliance_lag_mm  = NaN;
    end
end


%  SECTION 9: RUN — EARTH CONDITIONS
%  2 conditions (nominal, perturbed) x 2 controllers x 2 tests = 8 runs

fprintf('\n=== Running Earth Simulation ===\n');

[step_ref, step_vel] = make_trap_traj([0.3, 0.4], [0.7, 0.6], P.t_vec, ...
                                       P.vel_max, P.accel_max);

% Uncertainty model — applied only to perturbed condition
% Controller is never told; tests robustness to model mismatch
P.mass_uncertainty      = 0.15;
P.stiffness_uncertainty = 0.20;

nominal   = struct('mass_factor', 1.0, 'stiffness_factor', 1.0);
perturbed = struct('mass_factor', 1.0 + P.mass_uncertainty, ...
                   'stiffness_factor', 1.0 - P.stiffness_uncertainty);

conditions = struct(...
    'name',   {'Earth Nominal', 'Earth Perturbed'}, ...
    'g',      {GRAVITY.earth,   GRAVITY.earth},     ...
    'params', {nominal,         perturbed});

results = struct();

for c = 1:length(conditions)
    cond = conditions(c);
    fprintf('\n  Condition: %s (g=%.2f m/s^2)\n', cond.name, cond.g);

    for ctrl_idx = 1:2
        if ctrl_idx == 1
            cname = 'IT2';  fis_use = it2_fis_tuned;
        else
            cname = 'PID';  fis_use = [];
        end

        log1 = run_simulation(cname, fis_use, PID, P, cond.g, ...
                               step_ref, step_vel, false, cond.params);
        m1   = compute_metrics(log1, P);

        log2 = run_simulation(cname, fis_use, PID, P, cond.g, ...
                               step_ref, step_vel, true, cond.params);
        m2   = compute_metrics(log2, P);

        fname = sprintf('%s_%s', strrep(cond.name,' ','_'), cname);
        results.(fname).step         = log1;
        results.(fname).disturbed    = log2;
        results.(fname).metrics_step = m1;
        results.(fname).metrics_dist = m2;

        fprintf('    [%s] RMSE=%.2fmm | Settle=%.2fs | Recovery=%.2fs\n', ...
                cname, m1.rmse_pos, m1.settling_time, m2.dist_recovery_time);
    end
end


%  SECTION 9b: CONSOLE TABLE

fprintf('\n');
fprintf('%-24s | %-8s | %-9s | %-10s | %-10s | %-8s\n', ...
        'Condition + Controller', 'RMSE', 'MaxErr', 'Settle', 'Recovery', 'T_viol');
fprintf('%s\n', repmat('-', 1, 82));

cond_names = {'Earth_Nominal', 'Earth_Perturbed'};
ctrl_names = {'IT2', 'PID'};

for c = 1:length(cond_names)
    for ct = 1:length(ctrl_names)
        fname = sprintf('%s_%s', cond_names{c}, ctrl_names{ct});
        if isfield(results, fname)
            m1 = results.(fname).metrics_step;
            m2 = results.(fname).metrics_dist;
            fprintf('%-24s | %5.2fmm | %8.2fmm | %8.2fs | %8.2fs | %5d\n', ...
                    [strrep(cond_names{c},'_',' ') ' ' ctrl_names{ct}], ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, ...
                    m2.dist_recovery_time, m1.tension_violations);
        end
    end
end


%  SECTION 9c: FIGURES

figure('Name','CDPR Earth Results','Position',[50 50 1200 700]);

r_it2 = results.Earth_Nominal_IT2;
r_pid = results.Earth_Nominal_PID;

% 1. Step response
subplot(2,3,1);
plot(P.t_vec, sqrt(sum(r_it2.step.error.^2,2))*1000, 'b-', 'LineWidth',1.8); hold on;
plot(P.t_vec, sqrt(sum(r_pid.step.error.^2,2))*1000, 'r--','LineWidth',1.8);
yline(2,'k:','2% band'); grid on;
xlabel('Time (s)'); ylabel('Error (mm)');
title('Step Response: Tracking Error');
legend('IT2-FLS','PID','Location','NE');

% 2. Disturbance rejection
subplot(2,3,2);
plot(P.t_vec, sqrt(sum(results.Earth_Nominal_IT2.disturbed.error.^2,2))*1000, 'b-', 'LineWidth',1.8); hold on;
plot(P.t_vec, sqrt(sum(results.Earth_Nominal_PID.disturbed.error.^2,2))*1000, 'r--','LineWidth',1.8);
xline(P.disturbance_time,'k:','Disturbance');
xline(P.disturbance_time+P.disturbance_duration,'k:');
grid on; xlabel('Time (s)'); ylabel('Error (mm)');
title('Disturbance Rejection');
legend('IT2-FLS','PID');

% 3. Nominal vs perturbed (IT2)
subplot(2,3,3);
plot(P.t_vec, sqrt(sum(results.Earth_Nominal_IT2.step.error.^2,2))*1000, 'b-', 'LineWidth',1.8); hold on;
plot(P.t_vec, sqrt(sum(results.Earth_Perturbed_IT2.step.error.^2,2))*1000, 'g--','LineWidth',1.8);
grid on; xlabel('Time (s)'); ylabel('Error (mm)');
title('IT2: Nominal vs Perturbed');
legend('Nominal','Perturbed');

% 4. Nominal vs perturbed (PID)
subplot(2,3,4);
plot(P.t_vec, sqrt(sum(results.Earth_Nominal_PID.step.error.^2,2))*1000, 'r-', 'LineWidth',1.8); hold on;
plot(P.t_vec, sqrt(sum(results.Earth_Perturbed_PID.step.error.^2,2))*1000, 'm--','LineWidth',1.8);
grid on; xlabel('Time (s)'); ylabel('Error (mm)');
title('PID: Nominal vs Perturbed');
legend('Nominal','Perturbed');

% 5. FIS learning curve
subplot(2,3,5);
if ~isempty(cost_history)
    plot(1:length(cost_history), cost_history, 'b-', 'LineWidth',1.8);
    xlabel('Epoch'); ylabel('Cost J'); title('Neuro-Fuzzy Training'); grid on;
else
    text(0.5,0.5,'Fuzzy Toolbox unavailable','HorizontalAlignment','center','Units','normalized');
    title('Learning Curve (unavailable)');
end

% 6. RMSE bar chart
subplot(2,3,6);
rmse_vals = [results.Earth_Nominal_IT2.metrics_step.rmse_pos, ...
             results.Earth_Nominal_PID.metrics_step.rmse_pos, ...
             results.Earth_Perturbed_IT2.metrics_step.rmse_pos, ...
             results.Earth_Perturbed_PID.metrics_step.rmse_pos];
b = bar(rmse_vals, 'FaceColor','flat');
b.CData = [0.2 0.4 0.8; 0.9 0.3 0.3; 0.2 0.4 0.8; 0.9 0.3 0.3];
set(gca,'XTickLabel',{'Nom-IT2','Nom-PID','Pert-IT2','Pert-PID'},'XTickLabelRotation',25);
ylabel('RMSE (mm)'); title('RMSE by Condition'); grid on;

sgtitle('CDPR Earth Simulation: IT2-FLS vs PID');


%  SECTION 10: EXPORT

function export_results_table(results, P, filename_prefix)
    cond_names = {'Earth_Nominal', 'Earth_Perturbed'};
    ctrl_names = {'IT2', 'PID'};
    timestamp  = datestr(now, 'yyyy-mm-dd_HHMMSS');

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

    tex_file = [filename_prefix '_results_latex.tex'];
    fid = fopen(tex_file, 'w');
    fprintf(fid, '%% Auto-generated: %s\n', timestamp);
    fprintf(fid, '\\begin{table}[h]\n\\centering\n');
    fprintf(fid, '\\caption{CDPR Earth Simulation Performance Metrics}\n');
    fprintf(fid, '\\label{tab:earth_results}\n');
    fprintf(fid, '\\begin{tabular}{llrrrr}\n\\toprule\n');
    fprintf(fid, 'Condition & Controller & RMSE (mm) & Max Err (mm) & Settle (s) & Recovery (s) \\\\\n\\midrule\n');
    for c = 1:length(cond_names)
        for ct = 1:length(ctrl_names)
            fname = sprintf('%s_%s', cond_names{c}, ctrl_names{ct});
            if isfield(results, fname)
                m1 = results.(fname).metrics_step;
                m2 = results.(fname).metrics_dist;
                rt = m2.dist_recovery_time;
                if isnan(rt), rt_str = '---'; else, rt_str = sprintf('%.2f',rt); end
                fprintf(fid, '%s & %s & %.2f & %.2f & %.2f & %s \\\\\n', ...
                    strrep(cond_names{c},'_',' '), ctrl_names{ct}, ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, rt_str);
            end
        end
        fprintf(fid, '\\midrule\n');
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);

    fprintf('\n=== Export complete: %s, %s ===\n', csv_file, tex_file);
end

export_results_table(results, P, 'earth');

fprintf('\n=== Earth simulation complete. ===\n');
