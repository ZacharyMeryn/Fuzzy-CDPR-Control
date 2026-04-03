
% impactful asf
% Type-2 Neuro-Fuzzy vs PID Simulation Study — REDUCED GRAVITY CONDITIONS
% Cable-Driven Parallel Robot — 1m x 1m
%
% PURPOSE OF THIS FILE:
%   Run IT2-FLS vs PID across multiple reduced-gravity environments.
%   Gravity is just a number — swapping it is the entire mechanism for
%   simulating different bodies. No hardware changes needed.
%
%   Bodies simulated:
%     Moon    — 1.62  m/s^2  (1/6 Earth)  primary target for lunar construction
%     Mars    — 3.72  m/s^2  (0.38 Earth) next most relevant for exploration
%     Ceres   — 0.27  m/s^2  (0.03 Earth) smallest body, max sag stress test
%     Europa  — 1.32  m/s^2  (0.13 Earth) icy moon, interesting for future ops
%
%   Each body runs 2 controllers x 2 tests = 4 runs per body = 16 total.
%
% WHY THESE BODIES:
%   Moon and Mars are the obvious near-term targets. Ceres is the edge case —
%   almost zero gravity means cable sag is enormous and preload is critical.
%   Europa is a wildcard that shows the system behavior in a different low-g
%   regime. Together they give you a clean gravity sweep from 0.27 to 3.72 m/s^2
%   that makes for compelling plots.
%
% COMPANION FILE: ImpactfulFuzzballs_Earth_4_3_5pm.m
%   Runs the same controllers on Earth for the baseline comparison.
%
% BEFORE RUNNING ANYTHING, FILL IN MEASURED VALUES:
%   P.mass_EE          — weigh the end effector (kg)
%   P.cable_diameter   — measure with calipers (m)
%   P.cable_mass_per_m — from Dyneema spec sheet (kg/m)
%   P.cable_stiffness  — axial stiffness EA/L from spec sheet (N/m)
%   P.T_min / P.T_max  — from motor rated torque and cable break strength
%
% REQUIRED TOOLBOXES:
%   - Optimization Toolbox  (quadprog, lsqnonlin)
%   - Fuzzy Logic Toolbox   (R2019b+ for Type-2 support)

clear; clc; close all;


%  SECTION 0: PARAMETERS
%  Same robot parameters as the Earth file — only gravity changes per condition.

% Gravity values for each simulated body (m/s^2)
GRAVITY.earth  = 9.81;
GRAVITY.moon   = 1.62;
GRAVITY.mars   = 3.72;
GRAVITY.ceres  = 0.27;
GRAVITY.europa = 1.32;

% Workspace geometry
P.anchors = [0.0, 1.0;      % A1 top left
             1.0, 1.0;      % A2 top right
             1.0, 0.0;      % A3 bottom right
             0.0, 0.0];     % A4 bottom left
P.num_cables = 4;
P.ws_min = [0.05, 0.05];
P.ws_max = [0.95, 0.95];

% End-effector — placeholder, fill before running for real
P.mass_EE    = 0.5;
P.inertia_EE = 0.001;

% PTFE tube force
P.tube_force_per_meter = 0.3;

% Cable properties — placeholder, fill from spec sheet
P.cable_diameter  = 0.001;
P.cable_density   = 0.97;
P.cable_mass_per_m = P.cable_density * (pi*(P.cable_diameter/2)^2) * 1e6;
P.cable_stiffness = 50000;
P.cable_damping   = 5.0;

% Cable compliance bandwidth — needs professor check
P.cable_compliance_bandwidth = 20.0;   % Hz

% Tension limits
P.T_min = 2.0;
P.T_max = 50.0;

% Timing
P.dt      = 0.01;
P.t_total = 5.0;
P.t_vec   = 0:P.dt:P.t_total;
P.N       = length(P.t_vec);

% Disturbance
P.disturbance_time      = 2.5;
P.disturbance_magnitude = 1.5;
P.disturbance_duration  = 0.3;

% Trajectory limits
% >>> PROFESSOR CHECK: confirm before hardware runs
P.vel_max   = 0.20;
P.accel_max = 0.40;

fprintf('=== CDPR Low-Gravity Simulation Initialized ===\n');
fprintf('Bodies: Moon (%.2f) | Mars (%.2f) | Ceres (%.2f) | Europa (%.2f)\n\n', ...
        GRAVITY.moon, GRAVITY.mars, GRAVITY.ceres, GRAVITY.europa);
fprintf('Control rate: %.0f Hz | Duration: %.1f s | Steps: %d\n\n', ...
        1/P.dt, P.t_total, P.N);


%  SECTION 0b: TRAPEZOIDAL TRAJECTORY GENERATOR

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
            s = 0.5 * accel_max * t^2;   v = accel_max * t;
        elseif t <= t_end_flat
            s = 0.5 * accel_max * t_ramp^2 + v_peak * (t - t_end_ramp1);   v = v_peak;
        elseif t <= t_end_ramp2
            dt2 = t - t_end_flat;
            s   = 0.5 * accel_max * t_ramp^2 + v_peak * t_flat + v_peak * dt2 - 0.5 * accel_max * dt2^2;
            v   = v_peak - accel_max * dt2;
        else
            s = d;   v = 0;
        end
        s = min(s, d);
        pos_traj(k,:) = pos_start + s * dir;
        vel_traj(k,:) = v * dir;
    end
end


%  SECTION 1: INVERSE KINEMATICS WITH SAG

function [L_arc, sag_mid, L_chord] = ik_with_sag(pos, anchors, T_est, g, cable_mass_per_m)
    L_chord    = sqrt(sum((anchors - pos).^2, 2));
    w          = cable_mass_per_m * g;
    T_safe     = max(T_est, 0.1);
    sag_factor = (w .* L_chord).^2 ./ (24 .* T_safe.^2);
    L_arc      = L_chord .* (1 + sag_factor);
    sag_mid    = w .* L_chord.^2 ./ (8 .* T_safe);
end

% Sag comparison across all bodies at workspace center
T_demo   = ones(4,1) * 10;
pos_demo = [0.5, 0.5];
bodies   = {'moon','mars','ceres','europa'};
g_vals   = [GRAVITY.moon, GRAVITY.mars, GRAVITY.ceres, GRAVITY.europa];

fprintf('=== Cable Sag at Center (10N, Cable 1) ===\n');
fprintf('  %-10s  %-8s  %-12s\n', 'Body', 'g (m/s^2)', 'Sag (mm)');
for b = 1:length(bodies)
    [~, sag_b, ~] = ik_with_sag(pos_demo, P.anchors, T_demo, g_vals(b), P.cable_mass_per_m);
    fprintf('  %-10s  %-8.2f  %-12.4f\n', bodies{b}, g_vals(b), sag_b(1)*1e3);
end
fprintf('\n');


%  SECTION 2: FORWARD KINEMATICS

function [pos_est, valid, exit_flag] = forward_kinematics(L_measured, anchors, pos_guess)
    residual = @(p) sqrt(sum((anchors - p).^2, 2)) - L_measured;
    opts = optimoptions('lsqnonlin','Display','off', ...
        'FunctionTolerance',1e-10,'StepTolerance',1e-10,'MaxIterations',500);
    [pos_est, ~, ~, exit_flag] = lsqnonlin(residual, pos_guess, [0,0], [1,1], opts);
    valid = (exit_flag >= 1);
    if ~valid
        L_guess = sqrt(sum((anchors - pos_guess).^2, 2));
        J_ls    = (anchors - pos_guess) ./ L_guess;
        delta_p = (J_ls' * J_ls) \ (J_ls' * (L_measured - L_guess));
        pos_est = max([0,0], min([1,1], pos_guess + delta_p'));
        warning('FK:SolverFailed','lsqnonlin failed (flag=%d). Used WLS fallback.', exit_flag);
    end
end


%  SECTION 3: JACOBIAN AND TENSION DISTRIBUTION

function J = compute_jacobian(pos, anchors)
    L = sqrt(sum((anchors - pos).^2, 2));
    J = (anchors - pos) ./ L;
end

function [T, feasible] = tension_qp(pos, anchors, F_ext, T_min, T_max)
    J    = compute_jacobian(pos, anchors);
    n    = size(anchors, 1);
    opts = optimoptions('quadprog','Display','off');
    [T, ~, flag] = quadprog(eye(n), zeros(n,1), [], [], J', F_ext(:), ...
                            T_min*ones(n,1), T_max*ones(n,1), [], opts);
    feasible = (flag == 1);
    if ~feasible, T = (T_min + T_max)/2 * ones(n,1); end
end


%  SECTION 4: PID CONTROLLER
%  Per-body gain scheduling: lower gravity = lower cable tensions = reduced Kp.
%  Scaling is heuristic — proportional to sqrt(g/g_earth) as a first approximation
%  since cable restoring force scales roughly with tension, which scales with weight.
%
%  >>> PROFESSOR CHECK: gain scheduling law needs validation. sqrt(g/g_earth)
%  is physically motivated but not derived from the full closed-loop analysis.

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
    if F_mag > ctrl.output_lim, F_cmd = F_cmd * ctrl.output_lim / F_mag; end
    ctrl.e_prev = e;
end

function pid = make_pid_for_gravity(g, g_earth, dt)
    % Scale gains by sqrt(g/g_earth) — lower gravity = lower gains
    % Base gains tuned for Earth: Kp=18, Ki=0.8, Kd=3.0
    scale = sqrt(g / g_earth);
    Kp = 18 * scale;
    Ki = 0.8 * scale;
    Kd = 3.0 * scale;
    pid = pid_init([Kp, Kp], [Ki, Ki], [Kd, Kd], dt, 20.0);
end

% Build PID struct for each body
PID_earth  = make_pid_for_gravity(GRAVITY.earth,  GRAVITY.earth, P.dt);
PID_moon   = make_pid_for_gravity(GRAVITY.moon,   GRAVITY.earth, P.dt);
PID_mars   = make_pid_for_gravity(GRAVITY.mars,   GRAVITY.earth, P.dt);
PID_ceres  = make_pid_for_gravity(GRAVITY.ceres,  GRAVITY.earth, P.dt);
PID_europa = make_pid_for_gravity(GRAVITY.europa, GRAVITY.earth, P.dt);

fprintf('=== PID Gains by Body (Kp shown) ===\n');
fprintf('  Earth=%.2f | Moon=%.2f | Mars=%.2f | Ceres=%.2f | Europa=%.2f\n\n', ...
        PID_earth.Kp(1), PID_moon.Kp(1), PID_mars.Kp(1), ...
        PID_ceres.Kp(1), PID_europa.Kp(1));


%  SECTION 5: IT2 FUZZY LOGIC SYSTEM
%  Same architecture as Earth file. Uncertainty band is fixed at 15% for all bodies.
%  The FOU automatically provides more robustness at lower gravity where
%  the unmodeled dynamics (sag, compliance lag) are more pronounced.

function fis = build_it2_fls(sigma_uncertainty)
    fis = mamfis('Name','IT2_PositionController', ...
                 'AndMethod','min','OrMethod','max', ...
                 'ImplicationMethod','min','AggregationMethod','max', ...
                 'DefuzzificationMethod','centroid');
    fis = addInput(fis, [-0.4, 0.4], 'Name', 'pos_error');
    sig_e = 0.06;
    su    = sig_e * (1 + sigma_uncertainty);
    centers_e = [-0.25, -0.10, 0.0, 0.10, 0.25];
    mf_names  = {'NB','NS','ZE','PS','PB'};
    for i = 1:5
        fis = addMF(fis, 'pos_error', 'gaussmf', [su, centers_e(i)], 'Name', mf_names{i});
    end
    fis = addInput(fis, [-0.3, 0.3], 'Name', 'vel_error');
    sv_u = 0.05 * (1 + sigma_uncertainty);
    centers_v = [-0.20, -0.08, 0.0, 0.08, 0.20];
    for i = 1:5
        fis = addMF(fis, 'vel_error', 'gaussmf', [sv_u, centers_v(i)], 'Name', mf_names{i});
    end
    fis = addOutput(fis, [-20, 20], 'Name', 'force_cmd');
    out_centers = [-16, -8, 0, 8, 16];
    for i = 1:5
        fis = addMF(fis, 'force_cmd', 'trimf', ...
                    [out_centers(i)-4, out_centers(i), out_centers(i)+4], 'Name', mf_names{i});
    end
    rule_table = [1 1 1 2 2; 1 1 2 2 3; 1 2 3 4 5; 3 4 4 5 5; 4 4 5 5 5];
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


%  SECTION 6: NEURO-FUZZY TUNING
%  Tuned on Moon gravity since that's the primary target.
%  The Earth file tunes on Earth gravity.
%  Both share the same FIS architecture — only the MF centers differ after tuning.

function [fis_tuned, cost_history] = tune_fis_offline(fis, ref_trajectory, P, g, n_epochs, lr)
    if isempty(fis), fis_tuned = fis; cost_history = []; return; end
    fprintf('  Tuning FIS (%d epochs, g=%.2f)... ', n_epochs, g);
    cost_history = zeros(n_epochs, 1);
    fis_tuned    = fis;
    delta        = 1e-4;
    for epoch = 1:n_epochs
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
                new_c    = params(2) - lr * grad;
                in_range = fis_tuned.Inputs(inp).Range;
                new_c    = max(in_range(1)+0.01, min(in_range(2)-0.01, new_c));
                fis_tuned.Inputs(inp).MembershipFunctions(mf).Parameters(2) = new_c;
            end
        end
        if mod(epoch, 10) == 0, fprintf('.'); end
    end
    fprintf(' done. Cost: %.4f -> %.4f\n', cost_history(1), cost_history(end));
end

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

fprintf('\n=== Offline Neuro-Fuzzy Tuning (Moon gravity — primary target) ===\n');
if ~isempty(it2_fis)
    [it2_fis_tuned, cost_history] = tune_fis_offline(it2_fis, ref_traj_tune, ...
                                                       P, GRAVITY.moon, 30, 0.002);
else
    it2_fis_tuned = it2_fis;
    cost_history  = [];
end


%  SECTION 7: SIMULATION ENGINE

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
    error_mag = sqrt(sum(log.error.^2, 2)) * 1000;
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
        if isempty(recovery_idx)
            metrics.dist_recovery_time = NaN;
        else
            metrics.dist_recovery_time = recovery_idx * P.dt;
        end
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


%  SECTION 9: RUN — ALL REDUCED GRAVITY BODIES
%  4 bodies x 2 controllers x 2 tests = 16 runs
%  Each body uses nominal parameters only (no perturbed condition here —
%  that comparison lives in the Earth file).

fprintf('\n=== Running Low-Gravity Simulation Study ===\n');

[step_ref, step_vel] = make_trap_traj([0.3, 0.4], [0.7, 0.6], P.t_vec, ...
                                       P.vel_max, P.accel_max);

nominal = struct('mass_factor', 1.0, 'stiffness_factor', 1.0);

% Map body name to gravity value and PID struct
body_names = {'Moon', 'Mars', 'Ceres', 'Europa'};
body_g     = [GRAVITY.moon, GRAVITY.mars, GRAVITY.ceres, GRAVITY.europa];
body_pid   = {PID_moon, PID_mars, PID_ceres, PID_europa};

results = struct();

for b = 1:length(body_names)
    bname = body_names{b};
    g     = body_g(b);
    pid_b = body_pid{b};

    fprintf('\n  Body: %s (g=%.2f m/s^2)\n', bname, g);

    for ctrl_idx = 1:2
        if ctrl_idx == 1
            cname = 'IT2';  fis_use = it2_fis_tuned;
        else
            cname = 'PID';  fis_use = [];
        end

        log1 = run_simulation(cname, fis_use, pid_b, P, g, ...
                               step_ref, step_vel, false, nominal);
        m1   = compute_metrics(log1, P);

        log2 = run_simulation(cname, fis_use, pid_b, P, g, ...
                               step_ref, step_vel, true, nominal);
        m2   = compute_metrics(log2, P);

        fname = sprintf('%s_%s', bname, cname);
        results.(fname).step         = log1;
        results.(fname).disturbed    = log2;
        results.(fname).metrics_step = m1;
        results.(fname).metrics_dist = m2;
        results.(fname).g            = g;

        fprintf('    [%s] RMSE=%.2fmm | Settle=%.2fs | Recovery=%.2fs\n', ...
                cname, m1.rmse_pos, m1.settling_time, m2.dist_recovery_time);
    end
end


%  SECTION 9b: CONSOLE TABLE

fprintf('\n');
fprintf('%-20s | %-8s | %-9s | %-10s | %-10s | %-8s\n', ...
        'Body + Controller', 'RMSE', 'MaxErr', 'Settle', 'Recovery', 'T_viol');
fprintf('%s\n', repmat('-', 1, 80));

for b = 1:length(body_names)
    for ct = 1:length({'IT2','PID'})
        ctrl_names = {'IT2','PID'};
        fname = sprintf('%s_%s', body_names{b}, ctrl_names{ct});
        if isfield(results, fname)
            m1 = results.(fname).metrics_step;
            m2 = results.(fname).metrics_dist;
            fprintf('%-20s | %5.2fmm | %8.2fmm | %8.2fs | %8.2fs | %5d\n', ...
                    [body_names{b} ' ' ctrl_names{ct}], ...
                    m1.rmse_pos, m1.max_error, m1.settling_time, ...
                    m2.dist_recovery_time, m1.tension_violations);
        end
    end
end


%  SECTION 9c: FIGURES

figure('Name','CDPR Low-Gravity Results','Position',[50 50 1400 900]);

colors = {'b-','r--','g:','m-.'};   % one color per body

% 1. IT2 RMSE across all bodies — step response
subplot(2,3,1);
rmse_it2 = zeros(1,4);
rmse_pid = zeros(1,4);
for b = 1:4
    rmse_it2(b) = results.(sprintf('%s_IT2', body_names{b})).metrics_step.rmse_pos;
    rmse_pid(b) = results.(sprintf('%s_PID', body_names{b})).metrics_step.rmse_pos;
end
x = 1:4;
hold on;
bar(x-0.2, rmse_it2, 0.35, 'FaceColor',[0.2 0.4 0.8]);
bar(x+0.2, rmse_pid, 0.35, 'FaceColor',[0.9 0.3 0.3]);
set(gca,'XTick',1:4,'XTickLabel',body_names);
ylabel('RMSE (mm)'); title('Step Response RMSE by Body');
legend('IT2-FLS','PID'); grid on;

% 2. Settling time across bodies
subplot(2,3,2);
settle_it2 = zeros(1,4);
settle_pid = zeros(1,4);
for b = 1:4
    settle_it2(b) = results.(sprintf('%s_IT2', body_names{b})).metrics_step.settling_time;
    settle_pid(b) = results.(sprintf('%s_PID', body_names{b})).metrics_step.settling_time;
end
hold on;
bar(x-0.2, settle_it2, 0.35, 'FaceColor',[0.2 0.4 0.8]);
bar(x+0.2, settle_pid, 0.35, 'FaceColor',[0.9 0.3 0.3]);
set(gca,'XTick',1:4,'XTickLabel',body_names);
ylabel('Time (s)'); title('Settling Time by Body');
legend('IT2-FLS','PID'); grid on;

% 3. IT2 tracking error over time — all 4 bodies on one plot
subplot(2,3,3);
for b = 1:4
    fname = sprintf('%s_IT2', body_names{b});
    e = sqrt(sum(results.(fname).step.error.^2,2))*1000;
    plot(P.t_vec, e, colors{b}, 'LineWidth',1.8); hold on;
end
yline(2,'k:','2% band'); grid on;
xlabel('Time (s)'); ylabel('Error (mm)');
title('IT2 Tracking: All Bodies');
legend(body_names,'Location','NE');

% 4. PID tracking error over time — all 4 bodies
subplot(2,3,4);
for b = 1:4
    fname = sprintf('%s_PID', body_names{b});
    e = sqrt(sum(results.(fname).step.error.^2,2))*1000;
    plot(P.t_vec, e, colors{b}, 'LineWidth',1.8); hold on;
end
yline(2,'k:','2% band'); grid on;
xlabel('Time (s)'); ylabel('Error (mm)');
title('PID Tracking: All Bodies');
legend(body_names,'Location','NE');

% 5. Disturbance recovery — IT2 across bodies
subplot(2,3,5);
for b = 1:4
    fname = sprintf('%s_IT2', body_names{b});
    e = sqrt(sum(results.(fname).disturbed.error.^2,2))*1000;
    plot(P.t_vec, e, colors{b}, 'LineWidth',1.8); hold on;
end
xline(P.disturbance_time,'k:','Disturbance');
xline(P.disturbance_time+P.disturbance_duration,'k:');
grid on; xlabel('Time (s)'); ylabel('Error (mm)');
title('IT2 Disturbance Rejection: All Bodies');
legend(body_names,'Location','NE');

% 6. Cable sag at cable 1 — shows gravity effect clearly
subplot(2,3,6);
for b = 1:4
    fname = sprintf('%s_IT2', body_names{b});
    sag = results.(fname).step.sag_mid(:,1)*1000;
    plot(P.t_vec, sag, colors{b}, 'LineWidth',1.8); hold on;
end
grid on; xlabel('Time (s)'); ylabel('Sag at midspan (mm)');
title('Cable 1 Sag: All Bodies');
legend(body_names,'Location','NE');

sgtitle('CDPR Low-Gravity Study: IT2-FLS vs PID');


%  SECTION 10: EXPORT

function export_results_table(results, P, body_names, filename_prefix)
    ctrl_names = {'IT2', 'PID'};
    timestamp  = datestr(now, 'yyyy-mm-dd_HHMMSS');

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

    tex_file = [filename_prefix '_results_latex.tex'];
    fid = fopen(tex_file, 'w');
    fprintf(fid, '%% Auto-generated: %s\n', timestamp);
    fprintf(fid, '\\begin{table}[h]\n\\centering\n');
    fprintf(fid, '\\caption{CDPR Low-Gravity Simulation Performance Metrics}\n');
    fprintf(fid, '\\label{tab:lowgrav_results}\n');
    fprintf(fid, '\\begin{tabular}{llrrrrr}\n\\toprule\n');
    fprintf(fid, 'Body & Controller & g (m/s\\textsuperscript{2}) & RMSE (mm) & Max Err (mm) & Settle (s) & Recovery (s) \\\\\n\\midrule\n');
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
        fprintf(fid, '\\midrule\n');
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);

    fprintf('\n=== Export complete: %s, %s ===\n', csv_file, tex_file);
end

export_results_table(results, P, body_names, 'lowgrav');

fprintf('\n=== Low-gravity simulation complete. ===\n');
