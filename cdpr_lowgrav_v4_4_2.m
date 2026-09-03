% cdpr_lowgrav_v4_3.m
%
% v4.4.2 CHANGELOG -- removes a confound and a railed grid from v4.4.1.
%   (M28) NEW matched-Ki ablation. v4.4.1 gave the IT2 Ki = 800 and the T1
%         Ki = 1600, then showed the IT2 degrading less under mismatch
%         (+6.4% vs +10.6% on the Moon). A larger integral gain is
%         inherently more mismatch-sensitive, so that gap cannot be
%         attributed to the footprint until both run at the same Ki.
%         v4_matched_ki.csv reports exactly that comparison.
%   (M29) CFG.ki_grid extended to 6400: v4.4.1 railed the T1 and the SMC
%         at the 1600 ceiling, so neither of those designs is converged
%         and their numbers must not be published as optima.
%
% v4.4.1 CHANGELOG -- fixes three defects the v4.4 run exposed.
%   (M24) certification battery hold extended 2.5 -> 6.0 s whenever an
%         integrator is present. The tail p2p test was measuring an
%         unfinished integral transient and rejecting it as oscillation,
%         which is very likely why v4.4 selected Ki = 0 on three of four
%         bodies while Mars took Ki = 400. A certification window must
%         outlast the slowest closed-loop mode it is certifying.
%   (M25) tune_ki now prints the full scan (cost per gain, or BATT where
%         the battery rejected it). The v4.4 console reported only the
%         winner, making the body-to-body inconsistency unattributable.
%   (M26) the FOU sweep ran BOTH Ki = 0 and tuned Ki at each width. v4.4
%         retuned Ki per width, so the drop at delta = 0.30 coincided
%         exactly with Ki going 0 -> 400: footprint width and integral
%         gain moved together and the sweep no longer isolated either.
%   (M27) CFG.ki_grid extended to 1600; v4.4 railed at the 400 ceiling.
%
% v4.4 CHANGELOG -- matched integral action and a usable input scaling.
%   (M20) MATCHED INTEGRAL AUGMENTATION (CFG.integral_aug, default true).
%         In v4.3 the PID was the only controller with integral action.
%         Steady-state error against a constant bias is exactly what an
%         integrator removes, so that metric was decided structurally
%         before any simulation ran: fuzzy 0.53-0.58 mm and SMC 0.30 mm
%         against a PID essentially at zero. The fuzzy controllers and the
%         SMC now each carry an integral term with conditional integration
%         and anti-windup unwinding on saturation. Set false to reproduce
%         the v4.3 configuration exactly -- the with/without pair is a
%         clean ablation showing how much of the gap was architectural.
%   (M21) Ki selected in a second tuning stage with Ke/Ku (and lam/K/phi)
%         frozen, graded by the same dual-plant cost and required to clear
%         the same certification battery. Staged rather than joint so the
%         tuning cost stays additive instead of multiplying by the grid
%         size. Knee preference takes the SMALLEST adequate gain.
%   (M22) Ke_grid widened from [1 3 6 10 15] to [3 15 50 100 200 350].
%         The old grid was calibrated when tracking errors were 10-30 mm.
%         After the v4.3 limit-cycle fix errors fell to ~0.5 mm, which at
%         Ke = 3 normalizes to 0.0017 -- 0.4% of the input universe and
%         2.9% of one sigma. Membership activations there are
%         [0.001 0.337 1.000 0.363 0.002]: Zero saturated, NS and PS
%         firing almost symmetrically. The 25-rule base never left its
%         linear core, so the controller was a linear PD in all but name
%         and the FOU could contribute only its near-origin flattening
%         penalty. Reaching the NS/PS centres needs Ke of order 100-350.
%   (M23) the FOU sweep now mirrors the deployed configuration including
%         the integrator, or the delta = 0 and delta = 0.15 rows stop
%         reproducing Tables 1-2 and the ablation loses its consistency
%         check. v4_fou_sweep.csv gains a Ki column.
%
% v4.3.2 CHANGELOG -- output integrity and figure legibility.
%   (M15) root graphics defaults pinned explicitly. set(0,...) persists for
%         the MATLAB session, so a stale DefaultTextColor from an earlier
%         script silently recolors this one's labels with no trace in this
%         file. Most likely cause of the near-white callouts in the v4.3
%         configuration figure.
%   (M16) fig0 anchor/EE callouts drawn at FontSize 9, bold, explicit black.
%   (M17) NEW BLOCK 0 stale-output guard + end-of-run provenance manifest.
%         Outputs land in pwd; anything this run does not regenerate stays
%         behind looking current. A v4.2 FOU sweep CSV survived into a v4.3
%         folder this way and was one cross-check away from being cited.
%   (M18) fig8 uses the shared mm-based Springer sizing (174 mm x 70 mm)
%         instead of a pixel Position, matching fig0-fig7, and its bars now
%         take the shared CC greyscale palette. The v4.3.1 output rendered
%         in the MATLAB default colour order -- the only colour figure in
%         an otherwise greyscale set, and unusable in a print volume.
%   (M19) the effort verdict no longer fires on ratio alone. Two negligible
%         magnitudes can differ by orders of magnitude while meaning
%         nothing, so significance is gated on churn relative to the RMS
%         command, and the full-run ordering is printed because PID may
%         well be the smoothest commander of the four.
%
% v4.3.1 CHANGELOG -- adds control-effort / chattering instrumentation.
%   (M12) NEW BLOCK 9: control effort and chattering metrics computed from
%         the logs already held in R -- no re-simulation. Reports RMS
%         commanded force, command activity (mean |dF|/dt), normalized
%         total variation, command reversal rate, and realized-tension
%         rate, each over the full run AND over the station-keeping window
%         alone. Motivation: the four-way comparison scores tracking but
%         says nothing about what the switching costs. A boundary-layer
%         SMC buys its precision with command activity; without this
%         metric the comparison is silent on that trade.
%   (M13) lg.F is written at the CONTROL rate (100 Hz) but logged at the
%         PLANT rate (1 kHz), so it is held constant across P.decim
%         samples and 9 of every 10 increments are exactly zero. The
%         metrics are decimated by P.decim before differencing. NOTE: all
%         three rate metrics are duration-normalized and the reversal count
%         drops zero increments, so they are provably INVARIANT to the
%         logging rate (verified: ratio 0.998 undecimated vs decimated).
%         Decimation is therefore for clarity and speed, not correctness --
%         but it becomes load-bearing for any metric that is NOT
%         duration-normalized (e.g. a raw per-sample increment mean), so
%         the sequence is decimated once here and reused.
%   (M14) exports v4_control_effort.csv + fig8_control_effort.png.
%
% v4.3 CHANGELOG -- closes the remaining tuning blind spot found in the
% v4.2 run (IT2 nominal RMSE 1.46 mm was a mid-transit oscillation burst,
% MaxErr 6.7 mm at t ~ 2.3 s, invisible to the tail-only admissibility;
% the -70% "improvement" under perturbation was that burst switching off
% when the perturbed plant lost loop gain -- diagnostic, not robustness):
%   M6  CERTIFICATION BATTERY: oscillation admissibility is now also
%       checked by 2.5 s station-keeping holds (3 mm kick) at six poses
%       spanning the operating region, on BOTH plant parameterizations.
%       v4.2 registered IT2 onset at Ku = 1.0 from the hold-pose tail
%       (cap = 0.70) while transit poses oscillated below that: the
%       stability margin is pose-dependent through the structure matrix,
%       so the margin must be certified across poses, not at one pose.
%       Onset is bisection-refined (3 iters) before the 70% cap. The
%       battery gates PID and SMC selection identically (knee-preference
%       walk), and the FOU-sweep tuner inherits all of it.
%   M7  disturbance pulse moved 2.5 s -> 4.5 s: it previously fired
%       DURING the transit (motion ends at 2.74 s), superimposing the
%       tracking transient on the disturbance response in fig3.
%   M8  settling time measured against the FINAL pose. v4.2 measured
%       tracking error vs the MOVING reference, which never left the
%       2%-of-travel band (8.9 mm), so every Settle_s read 0.00. A 1 mm
%       absolute variant (Settle1mm_s) is also exported.
%   M9  recovery band = pre-disturbance error level + 0.5 mm. The old
%       2%-of-travel band exceeded every fuzzy/SMC disturbance peak
%       (recovery read 0.00 s) and barely trailed the PID peak (0.19 s).
%   M10 fig1(b) and Table 1 now report STEADY-STATE ERROR (settling is
%       trajectory-limited and near-identical across controllers after
%       M8). CSVs add RMSEhold_mm, Settle1mm_s, and the sweep's tuned Ku.
%   M11 legends moved off the bars (fig1, fig4); margin caps drawn on
%       fig6; deployed-FOU sanity note after the sweep (the v4.2 deployed
%       width 0.15 sat exactly on the sweep's marginal-stability spike).
% Runtime: ~45-75 min normal (battery adds ~50%), ~10-15 min fast_mode.
%
% v4.2 CHANGELOG -- implements the tuning procedure exactly as documented
% in manuscript Section 4.3 "Gain tuning", and adds the SMC baseline of
% Section 4.6. Fixes the RESIDUAL limit cycle of the v4.1 run:
%   M1  tuning trajectory now ends at the SAME hold pose as the
%       evaluation trajectory, (0.7, 0.6), via a different transit
%       (0.3, 0.9) -> (0.7, 0.6), motion complete at 3.0 s, 4 s of
%       scored station-keeping. v4.1 tuned at a pose that happened to be
%       barely oscillation-free while the evaluation pose was not
%       (stability margin is pose-dependent through the structure
%       matrix), so the oscillation penalty never fired.
%   M2  dual-plant cost: Eq. (11) summed over nominal AND perturbed
%       (mass x1.15, stiffness x0.85) plant parameterizations.
%   M3  explicit gain margin: within each Ke, output scales admissible
%       only up to 70% of the smallest scale exhibiting sustained
%       oscillation (tail p2p > 1.5 mm on either plant) or divergence;
%       knee rule then picks the smallest admissible scale within 5% of
%       minimum cost (>= 1.43x gain margin by construction).
%   M4  NEW fourth controller: boundary-layer SMC, F = K*sat(s/phi),
%       s = e_v + lam*e_p (switching only; gravity lives in the shared
%       allocation, no equivalent-control term). Identical cost,
%       admissibility, and least-aggressive knee.
%   M5  fig2 has 4 panels; tables/CSVs/figures carry 4 controllers.
% Runtime: ~25-40 min normal, ~6 min with CFG.fast_mode = true.
%
% v4.1 CHANGELOG (fixes the limit-cycle pathology found in the v4.0 run):
%   T1  tuning trajectory horizon 3 s -> 7 s. The tuning move takes 3.33 s,
%       so the v4.0 tuner never scored a single steady-state sample and
%       could not see the limit cycles it was creating.
%   T2  oscillation-penalized tuning cost: ISE + 25x tail ISE +
%       1000*(peak-to-peak tail error)^2. Sustained oscillation is now
%       far more expensive than a settled response.
%   T3  knee rule: select the SMALLEST gain within 5% of minimum cost
%       instead of the argmin (v4.0 selected the last stable point at the
%       base of the instability cliff -- see fig6 from the v4.0 run).
%   T4  input normalization gains Ke (error), Kv = Ke/2 (velocity): the
%       corrected plant produces 10-30 mm errors, the innermost ~5% of
%       the +/-0.4 m membership universe; Ke maps the operating range
%       onto the rule base (standard FLS Ke/Kv/Ku gain structure).
%   T5  PID selection knee: lowest bandwidth within 5% of minimum cost.
%   T6  disturbance recovery hold 0.5 s -> 2.0 s (a dip through the band
%       followed by divergence no longer counts as "recovered").
%   T7  NEW FOU width sweep on the Moon: fig7 + v4_fou_sweep.csv.
% =========================================================================
% IT2-FLS vs matched Type-1 FLS vs per-body-tuned PID
% Planar 4-cable CDPR, 1m x 1m, four reduced-gravity bodies
%
% v4 — full rebuild for the Springer revision. Every issue in
% ISSUES_MASTER_LIST.md items C1-C14 is fixed here:
%   C1  real IT2 engine: uncertain-sigma Gaussian FOU + Karnik-Mendel
%       type reduction (hand-rolled, matches revised Eq. (1)-(3) exactly).
%       Matched Type-1 baseline = same engine with FOU width 0.
%   C2  cable linear mass corrected (7.62e-4 kg/m, not 0.762 kg/m)
%   C3  two-rate sim: plant 1 kHz, controller+QP 100 Hz ZOH; exact
%       discrete first-order actuation lag (alpha ~ 0.118, no longer 1.0),
%       applied in the tension domain (tension-servoed winches); the
%       stiffness perturbation scales the loop bandwidth via tau = c/k
%   C4  allocation solves  A*tau = F_ctrl - w_hat  (nominal gravity
%       compensated in allocation; disturbances unknown to allocator)
%   C5  ONE shared force saturation (F_max) for every controller
%   C6  PID derivative acts on velocity error (same signal fuzzy sees)
%   C7  uniform velocity governor + inelastic workspace stops, applied
%       identically to all controllers
%   C8  perturbed condition (mass x1.15, stiffness x0.85) actually runs
%   C9  output-scale tuning by 1-D grid scan + local refine (well-posed)
%   C10 metrics: settling band defined (2% of travel), overshoot along
%       travel direction, recovery = band re-entry held >= 0.5 s
%   C12 base MATLAB (R2018b+) + Optimization Toolbox (quadprog) ONLY.
%       No Fuzzy Logic Toolbox, no Curve Fitting Toolbox.
%   C13 valid script layout: all local functions at end of file
%
% OUTPUTS (working directory):
%   v4_results_nominal.csv / v4_results_perturbed.csv
%   v4_table1_nominal.tex  / v4_table2_perturbed.tex
%   v4_tuning_summary.csv  (fuzzy scales + PID gains per body)
%   fig0_configuration.png     (NEW - reviewer item R8)
%   fig1_rmse_settling.png    (panel (b) = steady-state error since v4.3)
%   fig2_tracking_error.png    (3 panels: IT2 / T1 / PID, shared y)
%   fig3_disturbance_moon.png
%   fig4_robustness.png        (NEW centerpiece: nominal vs perturbed)
%   fig5_sag_lab_and_projection.png
%   fig6_scale_tuning.png      (replaces the old "learning curve")
%   fig7_fou_sweep.png + v4_fou_sweep.csv  (NEW: FOU width sweep, Moon)
%   fig8_control_effort.png + v4_control_effort.csv  (v4.3.1: effort/chatter)
%
% Runtime (v4.3): ~45-75 min normal, ~10-15 min with CFG.fast_mode = true.
% Self-checks run first and ABORT on physics errors. Do not comment
% them out; if one fails, something real is wrong.
% =========================================================================

clear; clc; close all;
t_wall = tic;

%% ===================== CONFIG =====================
CFG.fast_mode = false;     % true: dt_plant 2 ms + coarser tuning grids
CFG.fou_sweep = true;      % FOU width sweep on the Moon (fig7)
CFG.fou_deploy = 0.15;     % FOU width of the headline IT2 controller (see fig7 sanity note)
% v4.4 (M20): matched integral augmentation. In v4.3 the PID was the only
% controller with integral action, and steady-state error against a
% constant bias is exactly what integral action removes -- so that metric
% was structurally decided before any simulation ran. With this true the
% fuzzy controllers and the SMC each receive an integral term, selected by
% the same dual-plant cost and certification battery. false reproduces the
% v4.3 configuration exactly.
CFG.integral_aug = true;
CFG.ki_grid      = [0 25 50 100 200 400 800 1600 3200 6400];  % N/(m*s), after Ke/Ku
CFG.rng_seed  = 42;       % nothing stochastic in v4, kept for reproducibility
rng(CFG.rng_seed);

%% ===================== BLOCK 0: OUTPUT HYGIENE =====================
% v4.3.2 (M17): every export lands in pwd. If a previous version's outputs
% are already sitting there, any file this run does NOT regenerate stays
% behind looking current -- exactly how a v4.2 FOU sweep CSV survived into
% a v4.3 results folder and came within one cross-check of being cited.
% List them before anything is written, so the run cannot be confused with
% its predecessor.
stale = [dir('v4_*.csv'); dir('v4_*.tex'); dir('fig*.png')];
if ~isempty(stale)
    fprintf('\n[hygiene] %d pre-existing output file(s) in %s:\n', numel(stale), pwd);
    for i = 1:numel(stale)
        fprintf('   %-34s  %s\n', stale(i).name, datestr(stale(i).datenum, 'yyyy-mm-dd HH:MM'));
    end
    fprintf(['[hygiene] files NOT regenerated by this run will remain and will\n' ...
             '          look current. Clear the folder or use a fresh one.\n\n']);
end

%% ===================== BLOCK 1: PARAMETERS =====================
GRAV.moon = 1.62;  GRAV.mars = 3.72;  GRAV.ceres = 0.27;  GRAV.europa = 1.32;
body_names = {'Moon','Mars','Ceres','Europa'};
body_g     = [GRAV.moon, GRAV.mars, GRAV.ceres, GRAV.europa];
nb = numel(body_names);

P.anchors    = [0.0 1.0; 1.0 1.0; 1.0 0.0; 0.0 0.0];   % A1..A4
P.num_cables = 4;
P.ws_min = [0.05 0.05];   P.ws_max = [0.95 0.95];
P.mass_EE = 0.5;                       % kg (nominal; placeholder for hardware)
P.tube_force_per_meter = 0.3;          % N/m toward A1, unmodeled by controllers

% ---- C2 FIX: correct cable linear mass ----
P.cable_diameter   = 0.001;            % m
P.cable_density    = 970;              % kg/m^3 (Dyneema, 0.97 g/cm^3)
P.cable_mass_per_m = P.cable_density * pi*(P.cable_diameter/2)^2;  % 7.62e-4 kg/m

% Actuation compliance: winches are tension-controlled; realized tension
% follows commanded tension with a first-order lag whose bandwidth (20 Hz,
% Dyneema + servo, cf. Begey et al.) scales with cable stiffness since
% tau = c/k. The perturbed condition scales this bandwidth by stiff_f.
% (A literal elastic model T = T_cmd + k*dL with k = 50 kN/m saturates
% tension during any motion once the servo lag is resolvable - see
% README "expected behavior". Tension-domain lag is the coherent model.)
P.f_compliance    = 20.0;              % Hz, nominal tension-loop bandwidth

% ---- C3 FIX: two-rate timing ----
if CFG.fast_mode, P.dt_plant = 0.002; else, P.dt_plant = 0.001; end
P.dt_ctrl = 0.01;                                  % 100 Hz controller + QP
P.decim   = round(P.dt_ctrl / P.dt_plant);
P.alpha   = 1 - exp(-2*pi*P.f_compliance*P.dt_plant);   % exact discrete lag
P.t_total = 10.0;

P.T_min = 2.0;   P.T_max = 50.0;       % N
% v4.3 M7: pulse fires during station keeping. The eval transit ends at
% 2.74 s; the v4.2 pulse (2.5-2.8 s) fired DURING the move, superimposing
% the tracking transient on the disturbance response.
P.dist_t = 4.5;  P.dist_F = 1.5;  P.dist_dur = 0.3;
P.vel_max = 0.20;  P.accel_max = 0.40;

% ---- C5/C7 FIX: shared limits for every controller ----
P.F_max = 20.0;                        % N, shared saturation
P.K_gov = 40.0;                        % N*s/m velocity governor above vel_max

P.mass_unc  = 0.15;                    % perturbed condition (C8)
P.stiff_unc = 0.15;

fprintf('cdpr_lowgrav_v4_4_2 | plant %g Hz, control %g Hz, alpha=%.4f\n', ...
        1/P.dt_plant, 1/P.dt_ctrl, P.alpha);
fprintf('cable linear mass = %.3e kg/m (corrected)\n\n', P.cable_mass_per_m);

if exist('quadprog','file') ~= 2
    error('Optimization Toolbox (quadprog) is required.');
end

% evaluation trajectory: point-to-point, trapezoidal velocity (C14: not a "step")
t_eval = (0:P.dt_plant:P.t_total)';
[eval_ref, eval_vel] = make_trap_traj([0.3 0.4], [0.7 0.6], t_eval, P.vel_max, P.accel_max);

% tuning trajectory (v4.2 M1, manuscript Sec. 4.3): ends at the SAME hold
% pose as the evaluation trajectory via a different transit. d = 0.5 m,
% so motion completes at exactly 3.0 s; the 7 s horizon scores 4 s of
% station-keeping at the pose where station-keeping is actually graded.
t_tune = (0:P.dt_plant:7.0)';
[tune_ref, tune_vel] = make_trap_traj([0.3 0.9], [0.7 0.6], t_tune, P.vel_max, P.accel_max);

% v4.3 M6: certification-battery poses -- five samples along the
% evaluation transit plus the tuning-transit start. Stability margin is
% pose-dependent through the structure matrix, so gain admissibility is
% certified by station-keeping at each of these poses (see osc_battery).
frs = [0 0.25 0.5 0.75 1.0]';
P.cert_poses = [(1-frs)*[0.3 0.4] + frs*[0.7 0.6]; 0.3 0.9];

%% ===================== BLOCK 2: SELF-CHECKS (abort on failure) =====================
fprintf('self-checks: ');
run_selfchecks(P, GRAV.moon);
fprintf('all passed\n\n');

%% ===================== BLOCK 3: BUILD + TUNE CONTROLLERS =====================
fz_it2 = build_fls(CFG.fou_deploy);  % IT2: sigma_lo = (1-fou)*sigma, sigma_hi = (1+fou)*sigma
fz_t1  = build_fls(0.00);      % matched Type-1: FOU collapsed, same rules

if CFG.fast_mode
    grid_coarse = 0.2:0.8:6.0;  grid_fine_step = 0.10;  Ke_grid = [3 50 200];
else
    grid_coarse = 0.2:0.4:6.0;  grid_fine_step = 0.05;  Ke_grid = [3 15 50 100 200 350];
end
% If a tuned scale sits at the top of this range, the scan railed:
% widen the grid rather than publishing a boundary optimum.

sc_it2 = zeros(1,nb);  sc_t1 = zeros(1,nb);
ke_it2 = zeros(1,nb);  ke_t1 = zeros(1,nb);
scan_it2 = cell(1,nb); scan_t1 = cell(1,nb);
fprintf('fuzzy tuning (dual-plant cost, 70%% oscillation-onset cap, knee):\n');
for b = 1:nb
    fprintf('  %-7s ', body_names{b});
    [sc_it2(b), ke_it2(b), scan_it2{b}] = tune_fuzzy(fz_it2, P, body_g(b), tune_ref, tune_vel, grid_coarse, grid_fine_step, Ke_grid);
    fprintf('IT2: Ke=%g Ku=%.3f cap=%.2f ', ke_it2(b), sc_it2(b), scan_it2{b}.cap);
    [sc_t1(b), ke_t1(b), scan_t1{b}] = tune_fuzzy(fz_t1, P, body_g(b), tune_ref, tune_vel, grid_coarse, grid_fine_step, Ke_grid);
    fprintf('| T1: Ke=%g Ku=%.3f cap=%.2f\n', ke_t1(b), sc_t1(b), scan_t1{b}.cap);
end

if any([sc_it2 sc_t1] >= grid_coarse(end)-1e-9)
    warning('a fuzzy output scale railed at the top of the scan range - widen grid_coarse before publishing.');
end
fprintf('PID tuning (pole-placement grid, dual-plant cost + admissibility):\n');
pids = cell(1,nb);  pid_info = cell(1,nb);
for b = 1:nb
    [pids{b}, pid_info{b}] = tune_pid(P, body_g(b), tune_ref, tune_vel);
    fprintf('  %-7s f=%.1f Hz  zeta=%.1f  Ki-ratio=%.1f  (Kp=%.1f Ki=%.1f Kd=%.1f)\n', ...
        body_names{b}, pid_info{b}.f, pid_info{b}.z, pid_info{b}.kir, ...
        pids{b}.Kp, pids{b}.Ki, pids{b}.Kd);
end

fprintf('SMC tuning (lam x K x phi grid, dual-plant cost + admissibility):\n');
smcs = cell(1,nb);  smc_info = cell(1,nb);
for b = 1:nb
    [smcs{b}, smc_info{b}] = tune_smc(P, body_g(b), tune_ref, tune_vel);
    fprintf('  %-7s lam=%g 1/s  K=%g N  phi=%.2f m/s\n', ...
        body_names{b}, smcs{b}.lam, smcs{b}.K, smcs{b}.phi);
end

% ---- v4.4 (M21): integral-gain stage ----
% Two-stage rather than joint: Ke/Ku (and lam/K/phi) are selected first
% with the integrator off, then Ki is scanned with those frozen. A joint
% grid would multiply the tuning cost by numel(ki_grid); staging keeps it
% additive, and it mirrors standard practice of closing an inner loop
% before wrapping an outer integral around it. Ki is graded by the same
% dual-plant cost and must clear the same certification battery, so the
% "identical procedure for all four controllers" claim still holds.
ki_it2 = zeros(1,nb);  ki_t1 = zeros(1,nb);  ki_smc = zeros(1,nb);
if CFG.integral_aug
    fprintf('integral-gain stage (Ki scanned with Ke/Ku and lam/K/phi frozen):\n');
    for b = 1:nb
        base_it2 = struct('kind','fuzzy','fz',fz_it2,'scale',sc_it2(b),'Ke',ke_it2(b));
        base_t1  = struct('kind','fuzzy','fz',fz_t1, 'scale',sc_t1(b), 'Ke',ke_t1(b));
        fprintf('  %s\n', body_names{b});
        ki_it2(b) = tune_ki('IT2', base_it2, P, body_g(b), tune_ref, tune_vel, CFG.ki_grid, true);
        ki_t1(b)  = tune_ki('T1',  base_t1,  P, body_g(b), tune_ref, tune_vel, CFG.ki_grid, true);
        ki_smc(b) = tune_ki('SMC', smcs{b},  P, body_g(b), tune_ref, tune_vel, CFG.ki_grid, true);
        fprintf('    -> %-7s Ki: IT2=%g  T1=%g  SMC=%g  N/(m s)\n', ...
                body_names{b}, ki_it2(b), ki_t1(b), ki_smc(b));
    end
    if any([ki_it2 ki_t1 ki_smc] >= CFG.ki_grid(end)-1e-9)
        warning(['an integral gain railed at the top of CFG.ki_grid - widen it ' ...
                 'before publishing, the optimum is outside the scan.']);
    end
    for b = 1:nb, smcs{b}.Ki = ki_smc(b); end
else
    fprintf('integral augmentation DISABLED (CFG.integral_aug = false): v4.3 configuration.\n');
    for b = 1:nb, smcs{b}.Ki = 0; end
end

% tuning summary export
fid = fopen('v4_tuning_summary.csv','w');
fprintf(fid,'Body,IT2_Ke,IT2_Ku,IT2_Ki,T1_Ke,T1_Ku,T1_Ki,PID_f_Hz,PID_zeta,PID_Ki_ratio,Kp,Ki,Kd,SMC_lam,SMC_K,SMC_phi,SMC_Ki\n');
for b = 1:nb
    fprintf(fid,'%s,%g,%.4f,%g,%g,%.4f,%g,%.2f,%.2f,%.2f,%.4f,%.4f,%.4f,%g,%g,%.2f,%g\n', ...
        body_names{b}, ke_it2(b), sc_it2(b), ki_it2(b), ke_t1(b), sc_t1(b), ki_t1(b), ...
        pid_info{b}.f, pid_info{b}.z, pid_info{b}.kir, ...
        pids{b}.Kp, pids{b}.Ki, pids{b}.Kd, ...
        smcs{b}.lam, smcs{b}.K, smcs{b}.phi, ki_smc(b));
end
fclose(fid);

% ---- v4.4.2 (M28): matched-Ki ablation ----
% The v4.4.1 run selected Ki = 800 for the IT2 and Ki = 1600 for the T1,
% then reported the IT2 degrading less under mismatch (+6.4% vs +10.6% on
% the Moon). That comparison is CONFOUNDED: a larger integral gain is
% inherently more mismatch-sensitive, so the gap may be entirely the gain
% and not the footprint. Re-running both at a common Ki isolates the FOU,
% which is the whole point of keeping a matched Type-1 baseline.
if CFG.integral_aug
    ki_match = min([ki_it2 ki_t1]);
    fprintf(['matched-Ki ablation (both fuzzy controllers at Ki = %g, ' ...
             'isolating the FOU from the integral gain):\n'], ki_match);
    fid = fopen('v4_matched_ki.csv','w');
    fprintf(fid,'Body,Controller,Ki,RMSE_nom_mm,RMSE_pert_mm,dRMSE_pct,SS_mm\n');
    for b = 1:nb
        for cc = 1:2
            if cc == 1
                am = struct('kind','fuzzy','fz',fz_it2,'scale',sc_it2(b),'Ke',ke_it2(b),'Ki',ki_match);
                cn = 'IT2-FLS';  ct = 'IT2';
            else
                am = struct('kind','fuzzy','fz',fz_t1,'scale',sc_t1(b),'Ke',ke_t1(b),'Ki',ki_match);
                cn = 'T1-FLS';   ct = 'T1';
            end
            rm = zeros(1,2);
            for s = [1 3]
                lgm = run_sim(ct, am, P, body_g(b), eval_ref, eval_vel, scens{s});
                emm = sqrt(sum(lgm.err.^2,2))*1000;  emm = emm(isfinite(emm));
                if isempty(emm), rm(1 + (s==3)) = NaN; else, rm(1 + (s==3)) = sqrt(mean(emm.^2)); end
                if s == 1
                    ssm = mean(sqrt(sum(lgm.err(max(1,end-round(0.5/P.dt_plant)):end,:).^2,2)))*1000;
                end
            end
            fprintf(fid,'%s,%s,%g,%.4f,%.4f,%.2f,%.5f\n', body_names{b}, cn, ki_match, ...
                    rm(1), rm(2), (rm(2)-rm(1))/max(rm(1),eps)*100, ssm);
            fprintf('  %-7s %-8s RMSE %.4f -> %.4f mm (%+.1f%%)\n', ...
                    body_names{b}, cn, rm(1), rm(2), (rm(2)-rm(1))/max(rm(1),eps)*100);
        end
    end
    fclose(fid);
    fprintf('wrote v4_matched_ki.csv\n');
end

%% ===================== BLOCK 4: RUN MATRIX =====================
% 4 bodies x 4 controllers x 4 scenarios (C8: perturbation actually runs)
ctrl_list = {'IT2','T1','PID','SMC'};
clabels   = {'IT2-FLS','T1-FLS','PID','SMC'};
nc = numel(ctrl_list);

scens = { struct('name','nom_clean','inject',false,'mass_f',1.00,'stiff_f',1.00), ...
          struct('name','nom_dist', 'inject',true, 'mass_f',1.00,'stiff_f',1.00), ...
          struct('name','prt_clean','inject',false,'mass_f',1.15,'stiff_f',0.85), ...
          struct('name','prt_dist', 'inject',true, 'mass_f',1.15,'stiff_f',0.85) };

fprintf('\nrunning %d simulations...\n', nb*nc*numel(scens));
R = struct();  M = struct();
for b = 1:nb
    fprintf('  %s (g=%.2f)\n', body_names{b}, body_g(b));
    for c = 1:nc
        ck = ctrl_list{c};
        switch ck
            case 'IT2', asset = struct('kind','fuzzy','fz',fz_it2,'scale',sc_it2(b),'Ke',ke_it2(b),'Ki',ki_it2(b));
            case 'T1',  asset = struct('kind','fuzzy','fz',fz_t1, 'scale',sc_t1(b),'Ke',ke_t1(b),'Ki',ki_t1(b));
            case 'PID', asset = pids{b};
            case 'SMC', asset = smcs{b};
        end
        for s = 1:numel(scens)
            R.(body_names{b}).(ck).(scens{s}.name) = ...
                run_sim(ck, asset, P, body_g(b), eval_ref, eval_vel, scens{s});
        end
        M.(body_names{b}).(ck).nom = compute_metrics( ...
            R.(body_names{b}).(ck).nom_clean, R.(body_names{b}).(ck).nom_dist, P, eval_ref);
        M.(body_names{b}).(ck).prt = compute_metrics( ...
            R.(body_names{b}).(ck).prt_clean, R.(body_names{b}).(ck).prt_dist, P, eval_ref);
        mn = M.(body_names{b}).(ck).nom;  mp = M.(body_names{b}).(ck).prt;
        fprintf('    %-4s nom RMSE=%6.2f mm  ss=%5.2f mm | pert RMSE=%6.2f mm (%+5.1f%%)\n', ...
            ck, mn.rmse, mn.sse, mp.rmse, 100*(mp.rmse-mn.rmse)/max(mn.rmse,eps));
    end
end

%% ===================== BLOCK 5: CONSOLE TABLE + SANITY =====================
fprintf('\n%-14s %-6s | %8s %9s %8s %9s %9s | %8s %9s\n', ...
    'body','ctrl','RMSE(mm)','MaxE(mm)','SS(mm)','Overs(mm)','Recov(s)','pRMSE','dRMSE(%)');
fprintf('%s\n', repmat('-',1,104));
for b = 1:nb
    for c = 1:nc
        mn = M.(body_names{b}).(ctrl_list{c}).nom;
        mp = M.(body_names{b}).(ctrl_list{c}).prt;
        if isnan(mn.recovery), rs = '   --  '; else, rs = sprintf('%7.2f', mn.recovery); end
        fprintf('%-14s %-6s | %8.2f %9.2f %8.2f %9.2f %9s | %8.2f %8.1f%%\n', ...
            body_names{b}, ctrl_list{c}, mn.rmse, mn.maxe, mn.sse, ...
            mn.overshoot_mm, rs, mp.rmse, 100*(mp.rmse-mn.rmse)/max(mn.rmse,eps));
    end
end

% sanity warnings (do not hide these in the paper — investigate them)
for b = 1:nb
    for c = 1:nc
        mn = M.(body_names{b}).(ctrl_list{c}).nom;
        lg = R.(body_names{b}).(ctrl_list{c}).nom_clean;
        if mn.rmse > 100
            warning('%s %s: RMSE > 100 mm — investigate before publishing.', body_names{b}, ctrl_list{c});
        end
        if lg.qp_fail > 0
            warning('%s %s: %d QP infeasibilities.', body_names{b}, ctrl_list{c}, lg.qp_fail);
        end
    end
end

%% ===================== BLOCK 6: CSV + LATEX EXPORT =====================
export_csv('v4_results_nominal.csv',   M, body_names, body_g, ctrl_list, clabels, 'nom');
export_csv('v4_results_perturbed.csv', M, body_names, body_g, ctrl_list, clabels, 'prt');
export_tex_nominal('v4_table1_nominal.tex', M, body_names, body_g, ctrl_list, clabels);
export_tex_perturbed('v4_table2_perturbed.tex', M, body_names, ctrl_list, clabels);
fprintf('\nexported CSV + LaTeX tables.\n');

%% ===================== BLOCK 7: FIGURES =====================
mm2in = @(x) x/25.4;
CC = {[0 0 0],[0.45 0.45 0.45],[0.70 0.70 0.70],[0.20 0.20 0.20]};
LSty = {'-','--',':','-.'};
% v4.3.2 (M15): root graphics defaults set via set(0,...) PERSIST for the
% whole MATLAB session, so a value left behind by an earlier script can
% silently alter this one's output -- a faint DefaultTextColor, for
% instance, recolors every text() label without appearing anywhere in this
% file. Every default this script relies on is therefore pinned explicitly
% rather than inherited.
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',8, ...
      'DefaultTextFontName','Times New Roman','DefaultTextFontSize',8, ...
      'DefaultLegendFontSize',7,'DefaultLegendFontName','Times New Roman', ...
      'DefaultTextColor','k','DefaultAxesXColor','k','DefaultAxesYColor','k', ...
      'DefaultAxesColor','white','DefaultFigureColor','white', ...
      'DefaultLineLineWidth',1.0);

% ---- fig 0: configuration schematic (reviewer item R8) ----
fw = mm2in(84); fh = mm2in(80);
f0 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax = axes(f0); hold(ax,'on'); axis(ax,'equal');
ee0 = [0.55 0.50];
for i = 1:4
    plot(ax,[P.anchors(i,1) ee0(1)],[P.anchors(i,2) ee0(2)],'-','Color',[0.45 0.45 0.45],'LineWidth',1.0);
end
plot(ax,[P.ws_min(1) P.ws_max(1) P.ws_max(1) P.ws_min(1) P.ws_min(1)], ...
        [P.ws_min(2) P.ws_min(2) P.ws_max(2) P.ws_max(2) P.ws_min(2)], ...
        ':','Color',[0.6 0.6 0.6],'LineWidth',0.8);
plot(ax,P.anchors(:,1),P.anchors(:,2),'ks','MarkerSize',9,'MarkerFaceColor','k');
plot(ax,ee0(1),ee0(2),'ko','MarkerSize',9,'MarkerFaceColor','w','LineWidth',1.2);
labs = {'A_1','A_2','A_3','A_4'};
offs = [-0.10 0.05; 0.03 0.05; 0.03 -0.07; -0.10 -0.07];
for i = 1:4
    % v4.3.2 (M16): colour and weight stated explicitly. These labels
    % rendered near-white in the v4.3 output; production editors reject
    % schematics whose callouts do not survive greyscale printing.
    text(ax,P.anchors(i,1)+offs(i,1),P.anchors(i,2)+offs(i,2),labs{i}, ...
         'FontSize',9,'Color','k','FontWeight','bold');
end
text(ax,ee0(1)+0.04,ee0(2)-0.02,'EE','FontSize',9,'Color','k','FontWeight','bold');
% tube annotation along A1 cable
midp = 0.5*(P.anchors(1,:)+ee0);
text(ax,midp(1)-0.05,midp(2)+0.06,'feed tube','FontSize',8,'Color',[0.25 0.25 0.25]);
xlim(ax,[-0.18 1.18]); ylim(ax,[-0.15 1.15]);
springer_ax(ax,'Planar 4-cable CDPR configuration','x (m)','y (m)');
print(f0,'-dpng','-r600','fig0_configuration.png');

% ---- fig 1: RMSE + steady-state error grouped bars (v4.3 M10) ----
% Settling is trajectory-limited after the M8 fix (near-identical for all
% controllers); steady-state error is the differentiating station-keeping
% metric and feeds the KM-flattening discussion point.
RM = zeros(nb,nc); SS = zeros(nb,nc);
for b = 1:nb, for c = 1:nc
    RM(b,c) = M.(body_names{b}).(ctrl_list{c}).nom.rmse;
    SS(b,c) = M.(body_names{b}).(ctrl_list{c}).nom.sse;
end, end
fw = mm2in(174); fh = mm2in(70);
f1 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax1 = subplot(1,2,1); hb = bar(ax1,RM,'grouped');
for c = 1:nc, hb(c).FaceColor = CC{c}; hb(c).EdgeColor = 'none'; end
set(ax1,'XTick',1:nb,'XTickLabel',body_names);
lg = legend(ax1,clabels,'Location','northoutside','NumColumns',2,'Box','off'); fix_legend(lg);
springer_ax(ax1,'(a) Tracking RMSE','Planetary body','RMSE (mm)');
ax2 = subplot(1,2,2); hb = bar(ax2,SS,'grouped');
for c = 1:nc, hb(c).FaceColor = CC{c}; hb(c).EdgeColor = 'none'; end
set(ax2,'XTick',1:nb,'XTickLabel',body_names);
lg = legend(ax2,clabels,'Location','northoutside','NumColumns',2,'Box','off'); fix_legend(lg);
springer_ax(ax2,'(b) Steady-State Error','Planetary body','Error (mm)');
print(f1,'-dpng','-r600','fig1_rmse_settling.png');

% ---- fig 2: tracking error time series, one panel per controller, shared y ----
ymax = 0;
for b = 1:nb, for c = 1:nc
    e = err_mag_mm(R.(body_names{b}).(ctrl_list{c}).nom_clean);
    ymax = max(ymax, max(e(isfinite(e))));
end, end
y_sh = max(ymax*1.1, 1);
fw = mm2in(174); fh = mm2in(62);
f2 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ptag = {'(a) ','(b) ','(c) ','(d) '};
for c = 1:nc
    axc = subplot(1,4,c); hold(axc,'on');
    for b = 1:nb
        lg0 = R.(body_names{b}).(ctrl_list{c}).nom_clean;
        plot(axc, lg0.t, err_mag_mm(lg0), LSty{b}, 'Color', CC{b}, 'LineWidth', 1.1);
    end
    ylim(axc,[0 y_sh]);
    if c == 1, lg = legend(axc,body_names,'Location','northeast','Box','off'); fix_legend(lg); end
    springer_ax(axc,[ptag{c} clabels{c}],'Time (s)','Position error (mm)');
end
print(f2,'-dpng','-r600','fig2_tracking_error.png');

% ---- fig 3: disturbance rejection, Moon, all controllers ----
fw = mm2in(84); fh = mm2in(72);
f3 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax5 = axes(f3); hold(ax5,'on');
for c = 1:nc
    lg0 = R.Moon.(ctrl_list{c}).nom_dist;
    plot(ax5, lg0.t, err_mag_mm(lg0), LSty{c}, 'Color', CC{c}, 'LineWidth', 1.1);
end
xline(ax5,P.dist_t,':','Color',[0.5 0.5 0.5],'LineWidth',0.6,'HandleVisibility','off');
xline(ax5,P.dist_t+P.dist_dur,':','Color',[0.5 0.5 0.5],'LineWidth',0.6,'HandleVisibility','off');
lg = legend(ax5,clabels,'Location','northeast','Box','off'); fix_legend(lg);
springer_ax(ax5,'Disturbance Rejection (Moon)','Time (s)','Position error (mm)');
print(f3,'-dpng','-r600','fig3_disturbance_moon.png');

% ---- fig 4: robustness centerpiece — nominal vs perturbed RMSE ----
DG = zeros(nb,nc);
for b = 1:nb, for c = 1:nc
    mn = M.(body_names{b}).(ctrl_list{c}).nom.rmse;
    mp = M.(body_names{b}).(ctrl_list{c}).prt.rmse;
    DG(b,c) = 100*(mp-mn)/max(mn,eps);
end, end
fw = mm2in(84); fh = mm2in(72);
f4 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax6 = axes(f4); hb = bar(ax6,DG,'grouped');
for c = 1:nc, hb(c).FaceColor = CC{c}; hb(c).EdgeColor = 'none'; end
set(ax6,'XTick',1:nb,'XTickLabel',body_names);
lg = legend(ax6,clabels,'Location','northoutside','NumColumns',4,'Box','off'); fix_legend(lg);
springer_ax(ax6,'RMSE Degradation Under Parametric Mismatch','Planetary body','\Delta RMSE (%)');
print(f4,'-dpng','-r600','fig4_robustness.png');

% ---- fig 5: sag — (a) lab-scale simulated (µm), (b) construction-scale projection ----
fw = mm2in(174); fh = mm2in(66);
f5 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax7 = subplot(1,2,1); hold(ax7,'on');
for b = 1:nb
    lg0 = R.(body_names{b}).IT2.nom_clean;
    plot(ax7, lg0.t, movmean(lg0.sag(:,1),51)*1e6, LSty{b}, 'Color', CC{b}, 'LineWidth', 1.1);
end
lg = legend(ax7,body_names,'Location','southeast','Box','off'); fix_legend(lg);
springer_ax(ax7,'(a) Cable 1 Midspan Sag (laboratory scale)','Time (s)','Sag (\mum)');
ax8 = subplot(1,2,2); hold(ax8,'on');
span = linspace(5,50,200);
rho_big = 970*pi*0.003^2;   % 6 mm Dyneema, ~0.027 kg/m (illustrative scaled system)
T_big   = 200;              % N representative working tension
for b = 1:nb
    sagp = (rho_big*body_g(b)).*span.^2./(8*T_big);
    plot(ax8, span, sagp*1000, LSty{b}, 'Color', CC{b}, 'LineWidth', 1.1);
end
lg = legend(ax8,body_names,'Location','northwest','Box','off'); fix_legend(lg);
springer_ax(ax8,'(b) Projected Sag vs Span (6 mm cable, T = 200 N)','Span (m)','Sag (mm)');
print(f5,'-dpng','-r600','fig5_sag_lab_and_projection.png');

% ---- fig 6: output-scale tuning curves (replaces "learning curve") ----
fw = mm2in(84); fh = mm2in(65);
f6 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax9 = axes(f6); hold(ax9,'on');
s1 = scan_it2{1}; s2 = scan_t1{1};
plot(ax9, s1.s, s1.c, '-',  'Color', CC{1}, 'LineWidth', 1.1);
plot(ax9, s2.s, s2.c, '--', 'Color', CC{2}, 'LineWidth', 1.1);
i1 = find(abs(s1.s - s1.sel) < 1e-9, 1);  if isempty(i1), [~,i1] = min(s1.c); end
i2 = find(abs(s2.s - s2.sel) < 1e-9, 1);  if isempty(i2), [~,i2] = min(s2.c); end
plot(ax9, s1.s(i1), s1.c(i1), 'o', 'Color', CC{1}, 'MarkerFaceColor', CC{1}, 'MarkerSize', 4, 'HandleVisibility','off');
plot(ax9, s2.s(i2), s2.c(i2), 's', 'Color', CC{2}, 'MarkerFaceColor', CC{2}, 'MarkerSize', 4, 'HandleVisibility','off');
% v4.3 M11: show where the 70%-of-onset margin caps sit (battery-refined)
xline(ax9, s1.cap, ':', 'Color', CC{1}, 'LineWidth', 0.7, 'HandleVisibility','off');
xline(ax9, s2.cap, ':', 'Color', CC{2}, 'LineWidth', 0.7, 'HandleVisibility','off');
ylim(ax9, [0, 6*min([s1.c(:); s2.c(:)])]);
lg = legend(ax9,{'IT2-FLS','T1-FLS'},'Location','northwest','Box','off'); fix_legend(lg);
springer_ax(ax9,'Gain Tuning, Fine Stage (Moon)','Output scale K_u','Dual-plant cost \itJ');
print(f6,'-dpng','-r600','fig6_scale_tuning.png');

%% ===================== BLOCK 8: FOU WIDTH SWEEP (Moon) =====================
% v4.1 fix T7. Literature-motivated: reported IT2-vs-T1 deltas are known to
% depend strongly on FOU width, so the width itself is swept with all else
% fixed. Input gain held at the tuned IT2 Moon value; scale re-tuned per FOU.
if ~isfield(CFG,'fou_sweep') || CFG.fou_sweep
    fous = [0 0.15 0.30 0.45];
    fprintf('FOU sweep (Moon): ');
    fou_res = zeros(numel(fous), 6);   % [Ku, Ki, RMSE nom/pert at Ki=0, RMSE nom/pert at tuned Ki]
    for q = 1:numel(fous)
        fzq = build_fls(fous(q));
        sq = tune_ku(fzq, P, GRAV.moon, tune_ref, tune_vel, ...
                     grid_coarse, grid_fine_step, ke_it2(1));
        fou_res(q,1) = sq;
        % v4.4.1 (M26): v4.4 retuned Ki at each width, which CONFOUNDED the
        % sweep -- the drop at delta = 0.30 coincided with Ki jumping from 0
        % to 400, so FOU width and integral gain moved together and neither
        % could be credited. Both configurations are now run at every width:
        % Ki = 0 isolates the footprint, tuned Ki shows the deployed system.
        a0 = struct('kind','fuzzy','fz',fzq,'scale',sq,'Ke',ke_it2(1),'Ki',0);
        for s = [1 3]
            lgq = run_sim('IT2', a0, P, GRAV.moon, eval_ref, eval_vel, scens{s});
            emq = sqrt(sum(lgq.err.^2,2))*1000;
            fou_res(q, 3 + (s==3)) = sqrt(mean(emq(isfinite(emq)).^2));
        end
        aq = a0;
        if CFG.integral_aug
            aq.Ki = tune_ki('IT2', a0, P, GRAV.moon, tune_ref, tune_vel, CFG.ki_grid);
        end
        fou_res(q,2) = aq.Ki;
        for s = [1 3]
            lgq = run_sim('IT2', aq, P, GRAV.moon, eval_ref, eval_vel, scens{s});
            emq = sqrt(sum(lgq.err.^2,2))*1000;
            fou_res(q, 5 + (s==3)) = sqrt(mean(emq(isfinite(emq)).^2));
        end
        fprintf('.');
    end
    fprintf(' done\n');
    fid = fopen('v4_fou_sweep.csv','w');
    fprintf(fid,['FOU,Ku,Ki_tuned,RMSEnom_Ki0_mm,RMSEpert_Ki0_mm,' ...
                 'RMSEnom_Kituned_mm,RMSEpert_Kituned_mm\n']);
    for q = 1:numel(fous)
        fprintf(fid,'%.2f,%.4f,%g,%.4f,%.4f,%.4f,%.4f\n', fous(q), fou_res(q,1), ...
                fou_res(q,2), fou_res(q,3), fou_res(q,4), fou_res(q,5), fou_res(q,6));
    end
    fclose(fid);
    % v4.3 M11: deployed-FOU sanity note. In the v4.2 run the deployed
    % width (0.15) sat exactly on this sweep's marginal-stability spike.
    [~, qd] = min(abs(fous - CFG.fou_deploy));
    [bestR, qb] = min(fou_res(:,3));
    if fou_res(qd,3) > 1.10*bestR
        fprintf(['NOTE: deployed FOU %.2f -> nominal RMSE %.3f mm, but FOU %.2f -> %.3f mm.\n' ...
                 '      Consider CFG.fou_deploy = %.2f and a full rerun.\n'], ...
                 fous(qd), fou_res(qd,3), fous(qb), bestR, fous(qb));
    end
    fw = mm2in(84); fh = mm2in(65);
    f7 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
                'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
    ax10 = axes(f7); hold(ax10,'on');
    plot(ax10, fous, fou_res(:,3), '-o',  'Color', CC{1}, 'LineWidth', 1.1, ...
         'MarkerFaceColor', CC{1}, 'MarkerSize', 4);
    plot(ax10, fous, fou_res(:,4), '--s', 'Color', CC{2}, 'LineWidth', 1.1, ...
         'MarkerFaceColor', CC{2}, 'MarkerSize', 4);
    lg = legend(ax10,{'Nominal','Perturbed'},'Location','best','Box','off'); fix_legend(lg);
    springer_ax(ax10,'Effect of FOU Width (Moon, IT2-FLS)','FOU width','RMSE (mm)');
    print(f7,'-dpng','-r600','fig7_fou_sweep.png');
end

fprintf('figures saved at 600 DPI.\n');

%% ===================== BLOCK 9: CONTROL EFFORT / CHATTERING =====================
% Post-processing only: reads the logs already in R, runs no new simulation.
% Rationale: lg.F is the commanded planar control force, written once per
% control update (100 Hz) and held across the intervening plant steps
% (1 kHz), so 9 of every 10 plant-rate increments are exactly zero. Every
% metric below is evaluated on the decimated control-rate sequence. The
% three rate metrics happen to be invariant to this choice because they are
% normalized by duration and the reversal count discards zero increments,
% but decimating once up front keeps that independence explicit and makes
% the sequence safe to reuse for metrics that would NOT be invariant.
%
%   Frms      RMS commanded force magnitude (N)            -- gross effort
%   dFrate    mean |dF|/dt_ctrl over control steps (N/s)   -- command activity
%   TVn       total variation of F normalized by duration  -- chattering proxy
%   revs      sign reversals of dF per second, axis mean   -- direct chatter count
%   dTrate    mean |dT|/dt_ctrl on realized tensions (N/s) -- actuator wear proxy
%
% Each is reported over the full horizon AND over the station-keeping
% window alone (0.5 s after the reference stops moving to the end), since
% chattering while holding a pose is the operationally relevant case.

fprintf('\n=== CONTROL EFFORT / CHATTERING (nominal, clean) ===\n');

dref_ce = sqrt(sum(diff(eval_ref).^2, 2));
k_mv_ce = find(dref_ce > 1e-12, 1, 'last');   if isempty(k_mv_ce), k_mv_ce = 1; end
k_h_ce  = k_mv_ce + round(0.5 / P.dt_plant);  % start of the scored hold

EFF = struct();
fprintf('%-8s %-6s | %7s %9s %9s %8s %9s | %9s %9s\n', ...
    'body','ctrl','Frms(N)','dFrate','TVnorm','revs/s','dTrate','dFrate_h','revs/s_h');
fprintf('%s\n', repmat('-', 1, 88));

for b = 1:nb
    for c = 1:nc
        lg_ce = R.(body_names{b}).(ctrl_list{c}).nom_clean;

        idx_ce = 1:P.decim:size(lg_ce.F, 1);     % control-update instants
        Fc_ce  = lg_ce.F(idx_ce, :);
        Tc_ce  = lg_ce.T(idx_ce, :);
        ih_ce  = idx_ce >= k_h_ce;               % hold mask on the decimated grid

        e_full = eff_metrics(Fc_ce, Tc_ce, P.dt_ctrl);
        e_hold = eff_metrics(Fc_ce(ih_ce,:), Tc_ce(ih_ce,:), P.dt_ctrl);

        EFF.(body_names{b}).(ctrl_list{c}).full = e_full;
        EFF.(body_names{b}).(ctrl_list{c}).hold = e_hold;

        fprintf('%-8s %-6s | %7.3f %9.2f %9.2f %8.1f %9.2f | %9.2f %9.1f\n', ...
            body_names{b}, ctrl_list{c}, e_full.Frms, e_full.dFrate, ...
            e_full.TVn, e_full.revs, e_full.dTrate, e_hold.dFrate, e_hold.revs);
    end
end

export_effort('v4_control_effort.csv', EFF, body_names, ctrl_list, clabels);

% ---- figure: hold-window command activity and reversal rate ----
Ahold_ce = zeros(nb, nc);  Arev_ce = zeros(nb, nc);
for b = 1:nb
    for c = 1:nc
        Ahold_ce(b,c) = EFF.(body_names{b}).(ctrl_list{c}).hold.dFrate;
        Arev_ce(b,c)  = EFF.(body_names{b}).(ctrl_list{c}).hold.revs;
    end
end
fw = mm2in(174); fh = mm2in(70);
f8 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax11 = subplot(1,2,1); hb = bar(ax11, Ahold_ce, 'grouped');
for c = 1:nc, hb(c).FaceColor = CC{c}; hb(c).EdgeColor = 'none'; end
set(ax11,'XTick',1:nb,'XTickLabel',body_names);
lg = legend(ax11, clabels, 'Location','northoutside','NumColumns',2,'Box','off'); fix_legend(lg);
springer_ax(ax11,'(a) Command Activity (station keeping)','Planetary body','mean |\DeltaF|/\Deltat (N/s)');
ax12 = subplot(1,2,2); hb = bar(ax12, Arev_ce, 'grouped');
for c = 1:nc, hb(c).FaceColor = CC{c}; hb(c).EdgeColor = 'none'; end
set(ax12,'XTick',1:nb,'XTickLabel',body_names);
lg = legend(ax12, clabels, 'Location','northoutside','NumColumns',2,'Box','off'); fix_legend(lg);
springer_ax(ax12,'(b) Command Reversals (station keeping)','Planetary body','reversals per second');
print(f8,'-dpng','-r600','fig8_control_effort.png');
fprintf('wrote fig8_control_effort.png\n');

% ---- headline contrast for the manuscript (first body in the list) ----
b1 = body_names{1};
if isfield(EFF.(b1),'IT2') && isfield(EFF.(b1),'SMC') && isfield(EFF.(b1),'PID')
    fz_h = EFF.(b1).IT2.hold;  sm_h = EFF.(b1).SMC.hold;  pd_h = EFF.(b1).PID.hold;
    fprintf(['\n%s station keeping: IT2 %.2f N/s (%.1f rev/s) | ' ...
             'SMC %.2f N/s (%.1f rev/s) | PID %.2f N/s (%.1f rev/s)\n'], ...
        b1, fz_h.dFrate, fz_h.revs, sm_h.dFrate, sm_h.revs, pd_h.dFrate, pd_h.revs);
    fprintf('SMC/IT2 command-activity ratio during hold: %.2fx\n', ...
        sm_h.dFrate/max(fz_h.dFrate, eps));
    % v4.3.2 (M19): a ratio alone is not evidence. Two negligible numbers
    % can differ by orders of magnitude and mean nothing physically, so the
    % verdict also requires the activity to be a non-trivial fraction of
    % the commanded force itself before it is worth reporting as a cost.
    ratio_h = sm_h.dFrate / max(fz_h.dFrate, eps);
    % per-second command churn as a fraction of RMS command magnitude
    signif  = sm_h.dFrate / max(sm_h.Frms, eps);
    fprintf('  hold-window churn / RMS command: SMC %.2e per s, IT2 %.2e per s\n', ...
            signif, fz_h.dFrate/max(fz_h.Frms, eps));
    if ratio_h > 1.5 && signif > 0.05
        fprintf(['  -> SMC buys its precision with materially more command\n' ...
                 '     activity; report it as the cost axis.\n']);
    elseif ratio_h > 1.5
        fprintf(['  -> SMC is relatively noisier but ALL magnitudes are\n' ...
                 '     negligible (<5%% of RMS command per second). Report the\n' ...
                 '     ordering if useful, but do NOT frame it as a real cost.\n']);
    else
        fprintf(['  -> command activity is comparable; do NOT claim a chattering\n' ...
                 '     advantage for the fuzzy controllers.\n']);
    end
    % full-run ordering matters more than the hold alone, and PID may win it
    pd_f = EFF.(b1).PID.full;  fz_f = EFF.(b1).IT2.full;  sm_f = EFF.(b1).SMC.full;
    fprintf('  full-run dFrate: PID %.3f | IT2 %.3f | SMC %.3f N/s\n', ...
            pd_f.dFrate, fz_f.dFrate, sm_f.dFrate);
    if pd_f.dFrate < fz_f.dFrate
        fprintf(['  -> NOTE: PID is the smoothest commander overall. Do not\n' ...
                 '     present smoothness as a fuzzy advantage over PID.\n']);
    end
end

% ---- provenance manifest: what THIS run actually wrote ----
made = [dir('v4_*.csv'); dir('v4_*.tex'); dir('fig*.png')];
fprintf('\n[manifest] outputs in %s\n', pwd);
for i = 1:numel(made)
    age_min = (now - made(i).datenum)*24*60;
    if age_min < 120, tag = 'this run'; else, tag = 'STALE -- NOT from this run'; end
    fprintf('   %-34s  %s  %s\n', made(i).name, ...
            datestr(made(i).datenum,'yyyy-mm-dd HH:MM'), tag);
end

fprintf('\ntotal wall time: %.1f min\n', toc(t_wall)/60);
fprintf('done. numbers for the manuscript are in the v4_*.csv exports.\n');


% #########################################################################
% LOCAL FUNCTIONS (must remain at end of file)
% #########################################################################

function e = eff_metrics(F, T, dt)
    % Control-effort metrics on a CONTROL-RATE sequence (already decimated).
    % F: n x 2 commanded force; T: n x nc realized tensions; dt: control step.
    ok = all(isfinite(F), 2);
    F = F(ok, :);  T = T(ok, :);
    if size(F,1) < 3
        e = struct('Frms',NaN,'dFrate',NaN,'TVn',NaN,'revs',NaN,'dTrate',NaN);
        return;
    end
    dur  = (size(F,1)-1) * dt;
    Fmag = sqrt(sum(F.^2, 2));
    e.Frms = sqrt(mean(Fmag.^2));

    dF       = diff(F, 1, 1);
    dFmag    = sqrt(sum(dF.^2, 2));
    e.dFrate = mean(dFmag) / dt;              % N/s, mean command activity
    e.TVn    = sum(dFmag) / dur;              % total variation per second

    % sign reversals of the per-axis increment, averaged over axes. Zero
    % increments are dropped so that a held command does not register as a
    % reversal when it resumes in the opposite direction.
    r = 0;
    for a = 1:size(dF, 2)
        sg = sign(dF(:, a));
        sg = sg(sg ~= 0);
        if numel(sg) > 1, r = r + sum(diff(sg) ~= 0); end
    end
    e.revs = (r / size(dF, 2)) / dur;         % reversals per second per axis

    dT       = diff(T, 1, 1);
    e.dTrate = mean(sqrt(sum(dT.^2, 2))) / dt;  % N/s across the cable set
end

function export_effort(fname, EFF, body_names, ctrl_list, clabels)
    fid = fopen(fname, 'w');
    fprintf(fid, ['Body,Controller,Frms_N,dFrate_Nps,TVnorm_Nps,Reversals_ps,' ...
                  'dTrate_Nps,Frms_hold_N,dFrate_hold_Nps,Reversals_hold_ps\n']);
    for b = 1:numel(body_names)
        for c = 1:numel(ctrl_list)
            ef = EFF.(body_names{b}).(ctrl_list{c}).full;
            eh = EFF.(body_names{b}).(ctrl_list{c}).hold;
            fprintf(fid, '%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', ...
                body_names{b}, clabels{c}, ef.Frms, ef.dFrate, ef.TVn, ef.revs, ...
                ef.dTrate, eh.Frms, eh.dFrate, eh.revs);
        end
    end
    fclose(fid);
    fprintf('wrote %s\n', fname);
end

function run_selfchecks(P, g)
    % 1. discrete lag actually lags
    assert(P.alpha < 0.999, 'compliance lag saturated: reduce dt_plant or f_compliance');

    % 2. static hover: allocation with zero control must cancel gravity exactly
    pos = [0.5; 0.5];
    W_nom = [0; -P.mass_EE*g];
    [T, ok] = tension_alloc(pos, P.anchors, -W_nom, P.T_min, P.T_max);
    assert(ok, 'hover QP infeasible');
    Jm = jac(pos, P.anchors);
    res = Jm.'*T + W_nom;
    assert(norm(res) < 1e-6, 'hover force residual %.2e — allocation sign is wrong', norm(res));

    % 3. KM sanity: yl <= yr; FOU=0 collapses to exact weighted mean
    c = [-16 -8 0 8 16 -16 -8 0 8 16]';
    wl = rand(10,1); wh = wl + rand(10,1)*0.3;
    [yl, yr] = km_type_reduce(c, wl, wh);
    assert(yl <= yr + 1e-12, 'KM ordering violated');
    w = rand(10,1);
    [yl0, yr0] = km_type_reduce(c, w, w);
    ywm = sum(c.*w)/sum(w);
    assert(abs(yl0-ywm) < 1e-9 && abs(yr0-ywm) < 1e-9, 'KM does not collapse at FOU=0');

    % 4. IT2 engine with FOU=0 equals independent T1 weighted average
    fz0 = build_fls(0);
    for trial = 1:20
        e = (rand-0.5)*0.8;  v = (rand-0.5)*0.6;
        u_engine = fls_axis(fz0, e, v);
        mu_e = exp(-((e-fz0.e_c).^2)./(2*fz0.e_s^2)).';
        mu_v = exp(-((v-fz0.v_c).^2)./(2*fz0.v_s^2)).';
        w25 = min(mu_e(fz0.rule_i), mu_v(fz0.rule_j));
        u_ref = sum(fz0.rule_c.*w25)/max(sum(w25),1e-12);
        assert(abs(u_engine-u_ref) < 1e-9, 'IT2 engine mismatch at FOU=0');
    end

    % 5. corrected sag magnitude at lab scale is sub-mm
    [~, sag] = ik_sag([0.5;0.5], P.anchors, 10*ones(4,1), g, P.cable_mass_per_m);
    assert(max(sag) < 1e-3, 'lab-scale sag %.2e m — cable mass constant regressed?', max(sag));
end

% ------------------------- trajectory -------------------------
function [pos_traj, vel_traj] = make_trap_traj(p0, p1, t_vec, vmax, amax)
    N = numel(t_vec);  d = norm(p1-p0);  u = (p1-p0)/max(d,1e-9);
    vp = min(vmax, sqrt(amax*d));  tr = vp/amax;
    tf = (d - amax*tr^2)/vp;
    if tf < 0, tr = sqrt(d/amax); vp = amax*tr; tf = 0; end
    te1 = tr; te2 = tr+tf; te3 = 2*tr+tf;
    pos_traj = zeros(N,2); vel_traj = zeros(N,2);
    for k = 1:N
        t = t_vec(k);
        if t <= te1
            s = 0.5*amax*t^2;                              v = amax*t;
        elseif t <= te2
            s = 0.5*amax*tr^2 + vp*(t-te1);                v = vp;
        elseif t <= te3
            d2 = t-te2;
            s = 0.5*amax*tr^2 + vp*tf + vp*d2 - 0.5*amax*d2^2;
            v = vp - amax*d2;
        else
            s = d; v = 0;
        end
        s = min(s,d);
        pos_traj(k,:) = p0 + s*u;
        vel_traj(k,:) = v*u;
    end
end

% ------------------------- geometry -------------------------
function L = cab_len(pos, anchors)
    L = sqrt(sum((anchors - pos.').^2, 2));
end

function Jm = jac(pos, anchors)
    d = anchors - pos.';
    L = sqrt(sum(d.^2,2));
    Jm = d./L;                    % 4x2, rows are unit vectors toward anchors
end

function [L_arc, sag_mid] = ik_sag(pos, anchors, T, g, rho)
    Lc = cab_len(pos, anchors);
    w  = rho*g;
    Ts = max(T(:), 0.1);
    L_arc   = Lc.*(1 + (w.*Lc).^2 ./ (24*Ts.^2));
    sag_mid = w.*Lc.^2 ./ (8*Ts);
end

% ------------------------- allocation (C4 fix) -------------------------
function [T, ok] = tension_alloc(pos, anchors, F_des, Tmin, Tmax)
    persistent opts
    if isempty(opts), opts = optimoptions('quadprog','Display','off'); end
    n = size(anchors,1);
    if any(~isfinite(F_des(:)))
        T = 0.5*(Tmin+Tmax)*ones(n,1); ok = false; return;
    end
    Jm = jac(pos, anchors);
    [T,~,flag] = quadprog(eye(n), zeros(n,1), [], [], ...
                          Jm.', F_des(:), Tmin*ones(n,1), Tmax*ones(n,1), [], opts);
    ok = (flag == 1);
    if ~ok || isempty(T), T = 0.5*(Tmin+Tmax)*ones(n,1); end
end

% ------------------------- IT2 / T1 fuzzy engine (C1 fix) -------------------------
function fz = build_fls(fou)
    % uncertain-sigma Gaussian IT2 sets, revised Eq. (1):
    %   upper: exp(-(x-c)^2 / (2*sigma_hi^2)),  lower: same with sigma_lo,
    %   sigma_lo = sigma*(1-fou) <= sigma_hi = sigma*(1+fou).
    % fou = 0 gives the structurally matched Type-1 system.
    fz.fou = fou;
    fz.e_c = [-0.25 -0.10 0.0 0.10 0.25];  fz.e_s = 0.06;  fz.e_lim = 0.4;
    fz.v_c = [-0.20 -0.08 0.0 0.08 0.20];  fz.v_s = 0.05;  fz.v_lim = 0.3;
    out_c  = [-16 -8 0 8 16];
    rt = [1 1 1 2 2; 1 1 2 2 3; 1 2 3 4 5; 3 4 4 5 5; 4 4 5 5 5];
    [I, J2] = ndgrid(1:5, 1:5);            % column-major order matches rt(:)
    fz.rule_i = I(:);
    fz.rule_j = J2(:);
    fz.rule_c = reshape(out_c(rt(:)), [], 1);   % singleton consequents, 25x1
end

function u = fls_axis(fz, e, v)
    e = max(-fz.e_lim, min(fz.e_lim, e));
    v = max(-fz.v_lim, min(fz.v_lim, v));
    se_hi = fz.e_s*(1+fz.fou);  se_lo = max(fz.e_s*(1-fz.fou), 1e-9);
    sv_hi = fz.v_s*(1+fz.fou);  sv_lo = max(fz.v_s*(1-fz.fou), 1e-9);
    mu_e_hi = exp(-((e-fz.e_c).^2)./(2*se_hi^2)).';
    mu_e_lo = exp(-((e-fz.e_c).^2)./(2*se_lo^2)).';
    mu_v_hi = exp(-((v-fz.v_c).^2)./(2*sv_hi^2)).';
    mu_v_lo = exp(-((v-fz.v_c).^2)./(2*sv_lo^2)).';
    w_hi = min(mu_e_hi(fz.rule_i), mu_v_hi(fz.rule_j));
    w_lo = min(mu_e_lo(fz.rule_i), mu_v_lo(fz.rule_j));
    [yl, yr] = km_type_reduce(fz.rule_c, w_lo, w_hi);
    u = 0.5*(yl + yr);
end

function [yl, yr] = km_type_reduce(c, w_lo, w_hi)
    % Karnik-Mendel center-of-sets type reduction, singleton consequents
    c = c(:); w_lo = w_lo(:); w_hi = w_hi(:);
    if sum(w_hi) < 1e-12, yl = 0; yr = 0; return; end
    [c, ix] = sort(c);  w_lo = w_lo(ix);  w_hi = w_hi(ix);
    n = numel(c);
    % left endpoint: upper weights below switch point, lower above
    y = sum(c.*(w_lo+w_hi)) / sum(w_lo+w_hi);
    for it = 1:100
        k = find(c <= y, 1, 'last');
        if isempty(k), k = 1; end
        if k >= n, k = n-1; end
        w = [w_hi(1:k); w_lo(k+1:n)];
        yn = sum(c.*w)/max(sum(w),1e-12);
        if abs(yn - y) < 1e-10, y = yn; break; end
        y = yn;
    end
    yl = y;
    % right endpoint: lower weights below switch point, upper above
    y = sum(c.*(w_lo+w_hi)) / sum(w_lo+w_hi);
    for it = 1:100
        k = find(c <= y, 1, 'last');
        if isempty(k), k = 1; end
        if k >= n, k = n-1; end
        w = [w_lo(1:k); w_hi(k+1:n)];
        yn = sum(c.*w)/max(sum(w),1e-12);
        if abs(yn - y) < 1e-10, y = yn; break; end
        y = yn;
    end
    yr = y;
    if yl > yr, t = yl; yl = yr; yr = t; end
end

function F = fuzzy_ctrl(fz, e_pos, e_vel, scale, Ke)
    % v4.1 fix T4: input normalization (standard FLS Ke/Kv/Ku structure).
    % Corrected-plant errors are 10-30 mm, the innermost ~5% of the
    % +/-0.4 m design universe; Ke maps the operating range onto the
    % membership functions instead of forcing the output scale to
    % compensate. Kv = Ke/2 preserves the rule base's PD character
    % without saturating the velocity universe.
    Kv = max(1, Ke/2);
    e1 = Ke*e_pos;  v1 = Kv*e_vel;
    F = [fls_axis(fz, e1(1), v1(1)); fls_axis(fz, e1(2), v1(2))]*scale;
end

% ------------------------- PID (C5/C6 fixes) -------------------------
function C = pid_make(Kp, Ki, Kd, dt, F_max)
    C.kind = 'pid';
    C.Kp = Kp; C.Ki = Ki; C.Kd = Kd; C.dt = dt;
    C.int = [0;0];
    C.int_lim = 0.5*F_max/max(Ki,1e-9);     % anti-windup: I-force <= F_max/2
end

function [F, C] = pid_step(C, e_pos, e_vel)
    C.int = C.int + e_pos(:)*C.dt;
    C.int = max(-C.int_lim, min(C.int_lim, C.int));
    F = C.Kp*e_pos(:) + C.Ki*C.int + C.Kd*e_vel(:);
end

function [Cbest, info] = tune_pid(P, g, ref, refv)
    fgrid = [0.5 1.0 1.5 2.0 3.0];      % closed-loop bandwidth (Hz), << 20 Hz lag
    zgrid = [0.7 1.0 1.3];
    kgrid = [0.1 0.3 1.0];              % Ki as ratio of Kp
    % v4.2 (manuscript Sec. 4.3): identical dual-plant cost and
    % oscillation-admissibility filter; knee rule -- among admissible
    % designs within 5% of minimum cost, select the LOWEST bandwidth.
    P2P_TOL = 1.5e-3;
    n = 0; cand = struct('f',{},'z',{},'kir',{},'J',{},'p2p',{},'C',{});
    for a = 1:numel(fgrid)
        for z = 1:numel(zgrid)
            for q = 1:numel(kgrid)
                wc = 2*pi*fgrid(a);
                Kp = P.mass_EE*wc^2;
                Kd = 2*zgrid(z)*wc*P.mass_EE;
                Ki = kgrid(q)*Kp;
                C = pid_make(Kp, Ki, Kd, P.dt_ctrl, P.F_max);
                [J, pp] = dual_cost(C, 'PID', P, g, ref, refv);
                n = n + 1;
                cand(n) = struct('f',fgrid(a),'z',zgrid(z),'kir',kgrid(q), ...
                                 'J',J,'p2p',pp,'C',C);
            end
        end
    end
    adm = find([cand.p2p] <= P2P_TOL & [cand.J] < 1e9);
    if isempty(adm), error('tune_pid: no admissible design -- widen grids.'); end
    % v4.3 M6: identical certification battery gates the selection. Walk
    % the knee preference order (lowest bandwidth within 5% of minimum
    % cost, then every remaining tail-admissible design by cost) and take
    % the first design that also holds oscillation-free at every
    % certification pose on both plants.
    Jadm = [cand(adm).J];
    ok = adm(Jadm <= 1.05*min(Jadm));
    [~, ord] = sort([cand(ok).f] + 1e-6*[cand(ok).J]);
    rest = setdiff(adm, ok, 'stable');
    [~, o2] = sort([cand(rest).J]);
    walk = [ok(ord), rest(o2)];
    sel = 0;
    for w = walk
        [~, bad] = osc_battery('PID', cand(w).C, P, g, P2P_TOL);
        if ~bad, sel = w; break; end
    end
    if sel == 0, error('tune_pid: no design passes the certification battery -- widen grids.'); end
    Cbest = cand(sel).C;
    info = struct('f',cand(sel).f,'z',cand(sel).z,'kir',cand(sel).kir,'J',cand(sel).J);
end

% ------------------------- tuning (C9 fix) -------------------------
function [scale, Ke, scan] = tune_fuzzy(fz, P, g, ref, refv, coarse, fine_step, Ke_grid)
    % v4.2 M3 + v4.3 M6 (manuscript Sec. 4.3): per input gain Ke, scan Ku
    % ascending with the dual-plant cost until the oscillation onset
    % s_osc. Onset is now declared by EITHER criterion, on either plant:
    %   (i)  trajectory tail p2p > 1.5 mm or divergence (v4.2), or
    %   (ii) certification-battery failure: sustained oscillation or
    %        divergence during a 2.5 s station-keeping hold at ANY of the
    %        P.cert_poses (v4.3 -- the margin is pose-dependent through
    %        the structure matrix; the v4.2 tail-only criterion certified
    %        the hold pose while transit poses oscillated below the cap).
    % Onset is bisection-refined before the cap: Ku <= 0.7*s_osc.
    % Knee: smallest admissible Ku within 5% of the row minimum. Best
    % Ke = lowest knee-point cost. Fine pass repeats the tail criterion
    % inside the cap; the final pick is battery-verified.
    P2P_TOL = 1.5e-3;
    nK = numel(Ke_grid);
    rowJ   = nan(nK, numel(coarse));
    rowSel = nan(1, nK);  rowCost = inf(1, nK);  rowCap = inf(1, nK);
    for a = 1:nK
        s_osc = inf;
        for i = 1:numel(coarse)
            asset_i = fuzzy_asset(fz, coarse(i), Ke_grid(a));
            [J, pp] = dual_cost(asset_i, 'IT2', P, g, ref, refv);
            rowJ(a,i) = J;
            bad = (pp > P2P_TOL) || (J >= 1e9);
            if ~bad
                [~, bad] = osc_battery('IT2', asset_i, P, g, P2P_TOL);
            end
            if bad
                if i > 1, lo = coarse(i-1); else, lo = 0; end
                s_osc = refine_onset(@(u) fuzzy_asset(fz, u, Ke_grid(a)), ...
                                     'IT2', lo, coarse(i), P, g, P2P_TOL);
                break                              % early stop at onset
            end
        end
        cap = 0.7*s_osc;  rowCap(a) = cap;
        adm = find(coarse <= cap & isfinite(rowJ(a,:)));
        if isempty(adm), continue; end
        mn = min(rowJ(a,adm));
        ki = adm(find(rowJ(a,adm) <= 1.05*mn, 1, 'first'));
        rowSel(a) = coarse(ki);  rowCost(a) = rowJ(a,ki);
    end
    if all(isinf(rowCost))
        error('tune_fuzzy: no admissible gain found -- widen grids.');
    end
    [~, ba] = min(rowCost);
    Ke  = Ke_grid(ba);
    cap = rowCap(ba);
    if isinf(cap)
        warning('tune_fuzzy: no oscillation onset found in Ku range for Ke=%g; margin cap inactive.', Ke);
    end
    Ku0 = rowSel(ba);                              % battery-certified coarse knee
    fine = max(coarse(1), Ku0-0.35):fine_step:min(cap, Ku0+0.35);
    if isempty(fine), fine = Ku0; end
    fc = nan(size(fine));  fp = nan(size(fine));
    for i = 1:numel(fine)
        [fc(i), fp(i)] = dual_cost(fuzzy_asset(fz, fine(i), Ke), 'IT2', P, g, ref, refv);
    end
    adm = find(fp <= P2P_TOL & fc < 1e9);
    scale = Ku0;                                   % certified fallback
    while ~isempty(adm)
        mn = min(fc(adm));
        ii = adm(find(fc(adm) <= 1.05*mn, 1, 'first'));
        [~, bad] = osc_battery('IT2', fuzzy_asset(fz, fine(ii), Ke), P, g, P2P_TOL);
        if ~bad, scale = fine(ii); break; end
        adm(adm == ii) = [];                       % re-knee without it
    end
    scan.s = fine;  scan.c = fc;  scan.sel = scale;
    scan.coarse_s = coarse;  scan.coarse_J = rowJ;
    scan.Ke_grid = Ke_grid;  scan.cap = cap;  scan.Ke = Ke;
end

function Ki = tune_ki(ctype, base, P, g, ref, refv, ki_grid, verbose)
    % v4.4 (M21): select the integral gain with everything else frozen.
    % Graded by the same dual-plant cost as every other gain, and required
    % to clear the same station-keeping certification battery. Preference
    % is the knee: the SMALLEST gain whose cost is within 5% of the best,
    % because an oversized integrator buys steady-state accuracy with
    % overshoot and windup risk that the tail-weighted cost under-penalizes.
    P2P_TOL = 1.5e-3;
    nK = numel(ki_grid);
    J  = inf(1,nK);  ok = false(1,nK);
    for i = 1:nK
        a = base;  a.Ki = ki_grid(i);
        [jc, jp] = dual_cost(a, ctype, P, g, ref, refv);
        if ~isfinite(jc) || ~isfinite(jp), continue; end
        J(i) = jc + jp;
        [~, bad] = osc_battery(ctype, a, P, g, P2P_TOL);
        ok(i) = ~bad;
    end
    adm = find(ok & isfinite(J));
    % v4.4.1 (M25): report the full scan. In v4.4 the console showed only
    % the winner, so a selection of Ki = 400 on Mars and Ki = 0 on every
    % other body was impossible to attribute to cost or to the battery.
    if nargin >= 8 && verbose
        fprintf('    tune_ki(%-3s): ', ctype);
        for i = 1:nK
            if ~isfinite(J(i))
                fprintf('%g:div ', ki_grid(i));
            elseif ok(i)
                fprintf('%g:J=%.3e ', ki_grid(i), J(i));
            else
                fprintf('%g:BATT ', ki_grid(i));
            end
        end
        fprintf('\n');
    end
    if isempty(adm)
        warning('tune_ki(%s): no admissible integral gain; falling back to Ki = 0.', ctype);
        Ki = 0;  return;
    end
    knee = adm(J(adm) <= 1.05*min(J(adm)));
    Ki = min(ki_grid(knee));
end

function a = fuzzy_asset(fz, Ku, Ke)
    a = struct('kind','fuzzy','fz',fz,'scale',Ku,'Ke',Ke);
end

function [Ku, scan] = tune_ku(fz, P, g, ref, refv, coarse, fine_step, Ke)
    % Ku-only tuner for the FOU sweep: identical dual-plant cost, v4.3
    % combined admissibility (trajectory tail + certification battery),
    % bisection-refined onset, 70% cap, and knee, with Ke held fixed.
    P2P_TOL = 1.5e-3;
    rowJ = nan(size(coarse));  s_osc = inf;
    for i = 1:numel(coarse)
        asset_i = fuzzy_asset(fz, coarse(i), Ke);
        [J, pp] = dual_cost(asset_i, 'IT2', P, g, ref, refv);
        rowJ(i) = J;
        bad = (pp > P2P_TOL) || (J >= 1e9);
        if ~bad, [~, bad] = osc_battery('IT2', asset_i, P, g, P2P_TOL); end
        if bad
            if i > 1, lo = coarse(i-1); else, lo = 0; end
            s_osc = refine_onset(@(u) fuzzy_asset(fz, u, Ke), 'IT2', ...
                                 lo, coarse(i), P, g, P2P_TOL);
            break
        end
    end
    cap = 0.7*s_osc;
    adm = find(coarse <= cap & isfinite(rowJ));
    if isempty(adm), error('tune_ku: no admissible Ku -- widen grid.'); end
    mn = min(rowJ(adm));
    Ku0 = coarse(adm(find(rowJ(adm) <= 1.05*mn, 1, 'first')));
    fine = max(coarse(1), Ku0-0.35):fine_step:min(cap, Ku0+0.35);
    if isempty(fine), fine = Ku0; end
    fc = nan(size(fine));  fp = nan(size(fine));
    for i = 1:numel(fine)
        [fc(i), fp(i)] = dual_cost(fuzzy_asset(fz, fine(i), Ke), 'IT2', P, g, ref, refv);
    end
    adm = find(fp <= P2P_TOL & fc < 1e9);
    Ku = Ku0;                                      % certified fallback
    while ~isempty(adm)
        mn = min(fc(adm));
        ii = adm(find(fc(adm) <= 1.05*mn, 1, 'first'));
        [~, bad] = osc_battery('IT2', fuzzy_asset(fz, fine(ii), Ke), P, g, P2P_TOL);
        if ~bad, Ku = fine(ii); break; end
        adm(adm == ii) = [];
    end
    scan.s = fine;  scan.c = fc;  scan.sel = Ku;  scan.cap = cap;
end

% ---------------- v4.3 M6: certification battery ----------------
function [worst, failed] = osc_battery(ctype, asset, P, g, tol)
    % Station-keeping certification: 2.5 s holds at every pose in
    % P.cert_poses, on BOTH plant parameterizations, kicked by a 3 mm
    % initial offset (the hold reference steps from the kicked pose to
    % the true pose at the second sample; run_sim starts at ref(1,:)).
    % A design fails if the final-1.5 s peak-to-peak error exceeds tol at
    % any pose, or if any hold diverges (> 50 mm or non-finite).
    % Early-exits on the first failure.
    % v4.4.1 (M24): the certification window must outlast the slowest
    % closed-loop mode, or the tail measures an unfinished transient and
    % the p2p test reads slow integral settling as oscillation. In the
    % v4.4 run this rejected every Ki > 0 on three of four bodies. Holds
    % are extended whenever an integrator is present.
    has_int = (isfield(asset,'Ki') && isfinite(asset.Ki) && asset.Ki > 0) ...
              || strcmpi(ctype,'PID');
    if has_int, t_hold = 6.0; else, t_hold = 2.5; end
    t_tail = 1.5;
    N = round(t_hold/P.dt_plant) + 1;
    refv0 = zeros(N,2);
    kick  = 0.003*[1 1]/sqrt(2);
    scn = { struct('name','cert_nom','inject',false,'mass_f',1.00,'stiff_f',1.00), ...
            struct('name','cert_prt','inject',false,'mass_f',1.15,'stiff_f',0.85) };
    worst = 0;  failed = false;
    for q = 1:size(P.cert_poses,1)
        pose = P.cert_poses(q,:);
        r0 = min(max(pose + kick, P.ws_min + 0.002), P.ws_max - 0.002);
        ref0 = repmat(pose, N, 1);  ref0(1,:) = r0;
        for sc = 1:2
            lg = run_sim(ctype, asset, P, g, ref0, refv0, scn{sc});
            em = sqrt(sum(lg.err.^2,2));
            if any(~isfinite(em)) || max(em) > 0.05
                worst = inf;  failed = true;  return;
            end
            tail = em(end - round(t_tail/P.dt_plant):end);
            worst = max(worst, max(tail) - min(tail));
            if worst > tol, failed = true; return; end
        end
    end
end

function on = refine_onset(mk_asset, ctype, lo, hi, P, g, tol)
    % Bisect the oscillation onset between a battery-certified lo and a
    % failing hi (3 iterations, resolution ~ coarse_step/8). Returns the
    % smallest scale KNOWN to fail, so the 70% cap is conservative.
    for it = 1:3
        mid = 0.5*(lo + hi);
        [~, bad] = osc_battery(ctype, mk_asset(mid), P, g, tol);
        if bad, hi = mid; else, lo = mid; end
    end
    on = hi;
end

function [Cbest, info] = tune_smc(P, g, ref, refv)
    % v4.2 M4 (manuscript Sec. 4.6): boundary-layer SMC over the grid
    % lam x K x phi with the identical dual-plant cost and oscillation
    % admissibility. Knee: among admissible designs within 5% of minimum
    % cost, smallest switching gain K, then smallest lam, then largest
    % phi (least aggressive adequate design).
    P2P_TOL = 1.5e-3;
    lam_g = [2 4 8];  K_g = [2 5 10 20];  phi_g = [0.02 0.05 0.10];
    n = 0;  cand = struct('lam',{},'K',{},'phi',{},'J',{},'p2p',{},'C',{});
    for a = 1:numel(lam_g)
        for q = 1:numel(K_g)
            for w = 1:numel(phi_g)
                C = struct('kind','smc','lam',lam_g(a),'K',K_g(q),'phi',phi_g(w));
                [J, pp] = dual_cost(C, 'SMC', P, g, ref, refv);
                n = n + 1;
                cand(n) = struct('lam',lam_g(a),'K',K_g(q),'phi',phi_g(w), ...
                                 'J',J,'p2p',pp,'C',C);
            end
        end
    end
    adm = find([cand.p2p] <= P2P_TOL & [cand.J] < 1e9);
    if isempty(adm), error('tune_smc: no admissible design -- widen grids.'); end
    % v4.3 M6: identical certification battery gates the selection (same
    % walk as tune_pid, with the least-aggressive-adequate preference).
    Jadm = [cand(adm).J];
    ok = adm(Jadm <= 1.05*min(Jadm));
    score = [cand(ok).K]*1e4 + [cand(ok).lam]*1e2 - [cand(ok).phi];
    [~, ord] = sort(score);
    rest = setdiff(adm, ok, 'stable');
    [~, o2] = sort([cand(rest).J]);
    walk = [ok(ord), rest(o2)];
    sel = 0;
    for w = walk
        [~, bad] = osc_battery('SMC', cand(w).C, P, g, P2P_TOL);
        if ~bad, sel = w; break; end
    end
    if sel == 0, error('tune_smc: no design passes the certification battery -- widen grids.'); end
    Cbest = cand(sel).C;
    info = struct('lam',Cbest.lam,'K',Cbest.K,'phi',Cbest.phi,'J',cand(sel).J);
end

function [J, p2p] = traj_cost(asset, ctype, P, g, ref, refv, scen)
    % Eq. (11): ISE + 25x tail ISE (last 2 s) + 1000x (tail p2p)^2,
    % on the plant parameterization given by scen. Non-finite -> 1e9.
    lg = run_sim(ctype, asset, P, g, ref, refv, scen);
    if any(~isfinite(lg.err(:))), J = 1e9; p2p = inf; return; end
    em = sqrt(sum(lg.err.^2,2));
    J_ise  = sum(em.^2)*P.dt_plant;
    n_tail = min(numel(em)-1, round(2.0/P.dt_plant));
    tail   = em(end-n_tail:end);
    J_tail = sum(tail.^2)*P.dt_plant;
    p2p    = max(tail) - min(tail);
    J = J_ise + 25*J_tail + 1000*p2p^2;
    if ~isfinite(J), J = 1e9; p2p = inf; end
end

function [J, p2p] = dual_cost(asset, ctype, P, g, ref, refv)
    % v4.2 M2 (manuscript Sec. 4.3): Eq. (11) summed over the nominal and
    % perturbed plant parameterizations; p2p is the worst of the two and
    % feeds the oscillation-admissibility filter (M3).
    sn = struct('name','tune_nom','inject',false,'mass_f',1.00,'stiff_f',1.00);
    sp = struct('name','tune_prt','inject',false,'mass_f',1.15,'stiff_f',0.85);
    [Jn, pn] = traj_cost(asset, ctype, P, g, ref, refv, sn);
    [Jp, pp] = traj_cost(asset, ctype, P, g, ref, refv, sp);
    J = Jn + Jp;  p2p = max(pn, pp);
end

% ------------------------- simulation engine -------------------------
function lg = run_sim(ctype, asset, P, g, ref, refv, scen)
    % two-rate: plant at dt_plant, controller+allocation every P.decim steps
    N = size(ref,1);
    m_act = P.mass_EE * scen.mass_f;
    W_act = [0; -m_act*g];
    W_nom = [0; -P.mass_EE*g];              % allocator's belief (C4)
    % tension-loop lag; stiffness perturbation scales bandwidth (tau = c/k)
    f_act     = P.f_compliance * scen.stiff_f;
    alpha_act = 1 - exp(-2*pi*f_act*P.dt_plant);

    pid = [];
    if strcmpi(ctype,'PID'), pid = asset; pid.int = [0;0]; end
    % v4.4 (M20): matched integral augmentation for the non-PID controllers.
    % Ki = 0 reproduces v4.3 exactly. Integration is conditional: the state
    % does not accumulate while the command is saturated, the same
    % anti-windup discipline pid_step already applies.
    Ki_aug = 0;
    if isfield(asset,'Ki') && isfinite(asset.Ki), Ki_aug = asset.Ki; end
    e_int   = [0;0];
    int_cap = 0.5*P.F_max;          % |Ki*e_int| never exceeds half of saturation

    pos = ref(1,:).';  vel = [0;0];
    T_cmd   = tension_alloc(pos, P.anchors, -W_nom, P.T_min, P.T_max);
    T_state = T_cmd;                        % realized tension state
    F_ctrl = [0;0];

    lg.t = (0:N-1).' * P.dt_plant;
    lg.err = zeros(N,2);  lg.pos = zeros(N,2);  lg.vel = zeros(N,2);
    lg.F = zeros(N,2);    lg.T = zeros(N,P.num_cables);
    lg.sag = zeros(N,P.num_cables);
    lg.qp_fail = 0;  lg.sat_steps = 0;

    for k = 1:N
        t = lg.t(k);
        e_pos = ref(k,:).' - pos;
        e_vel = refv(k,:).' - vel;

        if mod(k-1, P.decim) == 0
            % ---- controller @ 100 Hz ----
            switch upper(ctype)
                case {'IT2','T1'}
                    F_ctrl = fuzzy_ctrl(asset.fz, e_pos, e_vel, asset.scale, asset.Ke);
                case 'PID'
                    [F_ctrl, pid] = pid_step(pid, e_pos, e_vel);
                case 'SMC'
                    % boundary-layer SMC (v4.2 M4, manuscript Sec. 4.6):
                    % switching term only -- gravity is compensated in the
                    % shared allocation, so no equivalent-control term.
                    s_smc  = e_vel + asset.lam*e_pos;
                    F_ctrl = asset.K * max(-1, min(1, s_smc/asset.phi));
                otherwise
                    F_ctrl = [0;0];
            end
            F_ctrl(~isfinite(F_ctrl)) = 0;
            % ---- matched integral term (fuzzy and SMC; PID has its own) ----
            if Ki_aug > 0 && ~strcmpi(ctype,'PID')
                e_int_try = e_int + e_pos*P.dt_ctrl;
                F_int     = Ki_aug*e_int_try;
                fi = norm(F_int);
                if fi <= int_cap                     % conditional integration
                    e_int = e_int_try;
                else
                    F_int = F_int*(int_cap/max(fi,eps));
                end
                F_ctrl = F_ctrl + F_int;
                F_ctrl(~isfinite(F_ctrl)) = 0;
            end
            % velocity governor (identical for all controllers, C7)
            sp = norm(vel);
            if sp > P.vel_max
                F_ctrl = F_ctrl - P.K_gov*(sp - P.vel_max)*(vel/sp);
            end
            % shared saturation (C5)
            fm = norm(F_ctrl);
            if fm > P.F_max
                F_ctrl = F_ctrl*(P.F_max/fm);
                lg.sat_steps = lg.sat_steps + 1;
                if Ki_aug > 0 && ~strcmpi(ctype,'PID')
                    e_int = e_int - e_pos*P.dt_ctrl;   % unwind while saturated
                end
            end
            % allocation with nominal gravity compensation (C4)
            [T_cmd, ok] = tension_alloc(pos, P.anchors, F_ctrl - W_nom, P.T_min, P.T_max);
            if ~ok, lg.qp_fail = lg.qp_fail + 1; end
        end

        % ---- plant @ 1 kHz: first-order tension-loop lag + bounds ----
        T_state = T_state + alpha_act*(T_cmd - T_state);
        T_act   = max(P.T_min, min(P.T_max, T_state));

        F_dist = [0;0];
        if scen.inject && t >= P.dist_t && t <= P.dist_t + P.dist_dur
            F_dist = [P.dist_F; 0];
        end
        d1 = P.anchors(1,:).' - pos;  nd1 = norm(d1);
        if nd1 > 1e-3
            F_dist = F_dist + P.tube_force_per_meter*nd1*(d1/nd1);
        end

        Jm = jac(pos, P.anchors);
        F_net = Jm.'*T_act + W_act + F_dist;
        acc = F_net/m_act;
        vel = vel + acc*P.dt_plant;                 % semi-implicit Euler
        pos_new = pos + vel*P.dt_plant;
        % workspace boundary = inelastic stop (C7)
        for ax = 1:2
            if pos_new(ax) < P.ws_min(ax)
                pos_new(ax) = P.ws_min(ax);
                if vel(ax) < 0, vel(ax) = 0; end
            elseif pos_new(ax) > P.ws_max(ax)
                pos_new(ax) = P.ws_max(ax);
                if vel(ax) > 0, vel(ax) = 0; end
            end
        end
        pos = pos_new;

        lg.err(k,:) = e_pos.';  lg.pos(k,:) = pos.';  lg.vel(k,:) = vel.';
        lg.F(k,:) = F_ctrl.';   lg.T(k,:) = T_act.';
        [~, sg] = ik_sag(pos, P.anchors, T_act, g, P.cable_mass_per_m);
        lg.sag(k,:) = sg.';

        if any(~isfinite(pos)) || any(~isfinite(vel))
            lg.err(k:end,:) = NaN;  break;
        end
    end
end

% ------------------------- metrics (C10 fix) -------------------------
function e = err_mag_mm(lg)
    e = sqrt(sum(lg.err.^2,2))*1000;
end

function m = compute_metrics(logClean, logDist, P, ref)
    em = err_mag_mm(logClean);
    ok = isfinite(em);
    if ~any(ok)
        m = struct('rmse',NaN,'rmse_hold',NaN,'maxe',NaN,'sse',NaN, ...
                   'settle',NaN,'settle1',NaN,'band',NaN,'overshoot_mm',NaN, ...
                   'recovery',NaN,'peak_dist_mm',NaN,'qp_fail',NaN,'stable',false);
        return;
    end
    em(~ok) = max(em(ok));
    m.rmse = sqrt(mean(em.^2));
    m.maxe = max(em);
    m.sse  = mean(em(max(1,end-round(0.5/P.dt_plant)):end));

    % v4.3: station-keeping RMSE (from 0.5 s after the reference stops
    % moving to the end of the horizon) -- separates the transit
    % transient from the hold, which total RMSE conflates.
    dref = sqrt(sum(diff(ref).^2,2));
    k_mv = find(dref > 1e-12, 1, 'last');  if isempty(k_mv), k_mv = 1; end
    k_h  = min(numel(em), k_mv + round(0.5/P.dt_plant));
    m.rmse_hold = sqrt(mean(em(k_h:end).^2));

    % v4.3 M8: settling measured against the FINAL pose. v4.2 measured
    % tracking error vs the MOVING reference, which never left the 2%
    % band (8.9 mm), so every settle read 0.00 s. Two variants: the 2%
    % band (trajectory-limited, near-identical across controllers) and a
    % 1 mm absolute band (differentiates the controllers).
    travel = norm(ref(end,:)-ref(1,:))*1000;       % mm
    m.band = 0.02*travel;                          % 2% of travel, mm
    d2f = sqrt(sum((logClean.pos - ref(end,:)).^2,2))*1000;
    m.settle  = settle_time(d2f, m.band, P.dt_plant, P.t_total);
    m.settle1 = settle_time(d2f, 1.0,    P.dt_plant, P.t_total);

    dirv = (ref(end,:)-ref(1,:))/max(norm(ref(end,:)-ref(1,:)),eps);
    proj = logClean.err*dirv.';                    % + behind target, - past it
    m.overshoot_mm = abs(max(0, -min(proj)*1000));   % abs() kills IEEE -0

    % v4.3 M9: recovery band relative to the PRE-disturbance error level
    % (baseline + 0.5 mm). The v4.2 band (2% of travel = 8.9 mm) was
    % wider than every fuzzy/SMC disturbance peak (recovery read 0.00 s)
    % and barely below the PID peak (0.19 s).
    ed = err_mag_mm(logDist);
    k0 = find(logDist.t >= P.dist_t, 1);
    m.recovery = NaN;  m.peak_dist_mm = NaN;
    hold_n = round(2.0/P.dt_plant);   % v4.1 fix T6
    if ~isempty(k0) && all(isfinite(ed))
        kb = max(1, k0 - round(0.5/P.dt_plant));
        rec_band = mean(ed(kb:k0)) + 0.5;          % mm
        for k = k0:(numel(ed)-hold_n)
            if all(ed(k:k+hold_n) < rec_band)
                m.recovery = (k-k0)*P.dt_plant;
                break;
            end
        end
        kend = min(k0 + round(2/P.dt_plant), numel(ed));
        m.peak_dist_mm = max(ed(k0:kend));
    end
    m.qp_fail = logClean.qp_fail;
    m.stable = all(isfinite(err_mag_mm(logClean))) && all(isfinite(err_mag_mm(logDist)));
end

function ts = settle_time(d, band, dt, t_total)
    % last time the distance-to-final-pose signal exceeds the band;
    % 0 if it never leaves the band, t_total if it never settles.
    idx = find(d > band);
    if isempty(idx)
        ts = 0;
    elseif idx(end) >= numel(d)
        ts = t_total;
    else
        ts = idx(end)*dt;
    end
end

% ------------------------- export -------------------------
function export_csv(fname, M, body_names, body_g, ctrl_list, clabels, which)
    fid = fopen(fname,'w');
    fprintf(fid,'Body,g_ms2,Controller,RMSE_mm,RMSEhold_mm,MaxErr_mm,SteadyState_mm,Settle_s,Settle1mm_s,Overshoot_mm,Recovery_s,PeakDistErr_mm,QPfail,Stable\n');
    for b = 1:numel(body_names)
        for c = 1:numel(ctrl_list)
            m = M.(body_names{b}).(ctrl_list{c}).(which);
            if isnan(m.recovery), rs = 'NaN'; else, rs = sprintf('%.3f', m.recovery); end
            fprintf(fid,'%s,%.2f,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%s,%.4f,%d,%d\n', ...
                body_names{b}, body_g(b), clabels{c}, m.rmse, m.rmse_hold, m.maxe, m.sse, ...
                m.settle, m.settle1, m.overshoot_mm, rs, m.peak_dist_mm, m.qp_fail, m.stable);
        end
    end
    fclose(fid);
end

function export_tex_nominal(fname, M, body_names, body_g, ctrl_list, clabels)
    fid = fopen(fname,'w');
    fprintf(fid,'%% auto-generated by cdpr_lowgrav_v4.m\n');
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Nominal-condition tracking performance. RMSE over the full 10 s trajectory; steady-state error: mean position error over the final 0.5 s of station keeping; recovery: time from disturbance onset until the error returns within 0.5 mm of its pre-disturbance level and remains there for 2.0 s.}\n');
    fprintf(fid,'\\label{tab:nominal}\n\\begin{tabular}{llrrrrr}\n\\toprule\n');
    fprintf(fid,'Body & Controller & $g$ (m/s$^2$) & RMSE (mm) & Max err (mm) & SS err (mm) & Recovery (s) \\\\\n\\midrule\n');
    for b = 1:numel(body_names)
        for c = 1:numel(ctrl_list)
            m = M.(body_names{b}).(ctrl_list{c}).nom;
            if isnan(m.recovery), rs = '---'; else, rs = sprintf('%.2f', m.recovery); end
            fprintf(fid,'%s & %s & %.2f & %.2f & %.2f & %.2f & %s \\\\\n', ...
                body_names{b}, clabels{c}, body_g(b), m.rmse, m.maxe, m.sse, rs);
        end
        if b < numel(body_names), fprintf(fid,'\\midrule\n'); end
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
end

function export_tex_perturbed(fname, M, body_names, ctrl_list, clabels)
    fid = fopen(fname,'w');
    fprintf(fid,'%% auto-generated by cdpr_lowgrav_v4.m\n');
    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,'\\caption{Robustness under parametric mismatch (mass $\\times$1.15, stiffness $\\times$0.85, unknown to all controllers).}\n');
    fprintf(fid,'\\label{tab:perturbed}\n\\begin{tabular}{llrrr}\n\\toprule\n');
    fprintf(fid,'Body & Controller & RMSE$_\\mathrm{pert}$ (mm) & $\\Delta$RMSE (\\%%) & Recovery$_\\mathrm{pert}$ (s) \\\\\n\\midrule\n');
    for b = 1:numel(body_names)
        for c = 1:numel(ctrl_list)
            mn = M.(body_names{b}).(ctrl_list{c}).nom;
            mp = M.(body_names{b}).(ctrl_list{c}).prt;
            dg = 100*(mp.rmse - mn.rmse)/max(mn.rmse, eps);
            if isnan(mp.recovery), rs = '---'; else, rs = sprintf('%.2f', mp.recovery); end
            fprintf(fid,'%s & %s & %.2f & %+.1f & %s \\\\\n', ...
                body_names{b}, clabels{c}, mp.rmse, dg, rs);
        end
        if b < numel(body_names), fprintf(fid,'\\midrule\n'); end
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
end

% ------------------------- figure helpers -------------------------
function springer_ax(ax, ttl, xl, yl)
    set(ax,'Box','off','TickDir','in','LineWidth',0.5,'Color','white', ...
           'XColor','k','YColor','k','GridColor',[0.85 0.85 0.85], ...
           'GridAlpha',1.0,'FontName','Times New Roman','FontSize',8);
    grid(ax,'on');
    if ~isempty(ttl), title(ax,ttl,'FontSize',9,'FontWeight','bold','FontName','Times New Roman','Color','k'); end
    if ~isempty(xl), xlabel(ax,xl,'FontSize',8,'FontName','Times New Roman','Color','k'); end
    if ~isempty(yl), ylabel(ax,yl,'FontSize',8,'FontName','Times New Roman','Color','k'); end
end

function fix_legend(lg)
    set(lg,'TextColor','k','EdgeColor',[0.5 0.5 0.5]);
end
