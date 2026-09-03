%% ========================================================================
%  cdpr_dynamic_disturbance.m  --  How do the four controllers hold up
%  against a TIME-VARYING disturbance they were never tuned against?
%
%  MOTIVATION
%  The v4.4.2 study exercises one impulsive disturbance and one constant
%  bias. Both are quasi-static and MATCHED (they enter through the same
%  channel as the control input), which is the ideal operating case for a
%  boundary-layer sliding mode controller: inside the layer, sat(s/phi)
%  acts as a proportional term of effective gain K/phi = 100 N per unit
%  surface, and that stiffness is what defeats a constant bias. The SMC
%  therefore wins that study on its home ground.
%
%  This file adds disturbance classes the plant does not currently
%  exercise, applied in EVALUATION ONLY with the gains already selected
%  and published. Nothing is retuned. That makes this a test of robustness
%  to an UNMODELLED disturbance class, not a competition to see which
%  controller adapts best -- and it means no controller can be accused of
%  having tuned its way to the result.
%
%  PREDICTION, STATED BEFORE THE RUN (record it, report it either way):
%    The same high effective gain that wins the static case amplifies
%    whatever it is fed. Under 'stribeck' the boundary layer sits exactly
%    where stick-slip is worst -- low velocity -- so the SMC is expected to
%    degrade MOST, the fuzzy controllers LEAST (bounded smooth surface, no
%    gain amplification), with the PID between them. If the ordering comes
%    out otherwise, that is the result and it gets reported.
%
%  DISTURBANCE MODELS (DIST.modes)
%   'stribeck'  Dust-contaminated winch friction, lumped to the EE. Coulomb
%               + Stribeck dip + viscous, opposing motion, discontinuous at
%               zero velocity. Motivated directly by Sect. 3.2: abrasive,
%               electrostatically adhesive regolith introducing systematic
%               friction and backlash. THE DEFENSIBLE DEFAULT.
%   'extrusion' Deposition reaction force: periodic pulsation from the
%               extruder plus a stochastic component. Physically the actual
%               application (MMPACT-class construction).
%   'gust'      Band-limited stochastic loading. ONLY defensible on Mars.
%               The Moon, Ceres and Europa have no atmosphere -- do not
%               report this for those bodies. The script refuses to run it
%               on airless bodies unless you override DIST.allow_airless.
%
%  PREREQUISITE: run cdpr_lowgrav_v4_4_2.m to completion first, then run
%  this in the same session. It reuses that workspace and auto-extracts the
%  v4.4.2 local functions, so the two files cannot drift apart.
%
%  OUTPUTS: v4_dyn_disturbance.csv, fig11_dyn_disturbance.png
%% ========================================================================

clear DIST;
DIST.modes       = {'stribeck','extrusion','gust'};  % all three in one run;
                                   % trim the list to run a subset
DIST.levels      = [0 0.5 1.0 2.0];% severity multipliers (0 = baseline)
DIST.n_seeds     = 6;              % replicates (stochastic modes only)
DIST.src         = 'cdpr_lowgrav_v4_4_2.m';
DIST.libdir      = fullfile(pwd, 'v442_lib');
DIST.allow_airless = false;        % guard against reporting 'gust' off-Mars

% --- Stribeck friction parameters (lumped to the EE) ---
% F_fric = -(Fc + (Fs-Fc)*exp(-(|v|/vs)^2) + sigma*|v|) * sign(v)
% Fc/Fs chosen as a few percent of the ~0.24 N RMS command, so the
% disturbance is a perturbation rather than a plant redesign.
DIST.Fc = 0.010;      % N, Coulomb level at severity 1.0
DIST.Fs = 0.018;      % N, static/breakaway level
DIST.vs = 0.004;      % m/s, Stribeck velocity
DIST.sig= 0.020;      % N/(m/s), viscous term

% --- Extrusion reaction parameters ---
DIST.ex_f    = 7.0;   % Hz, extruder pulsation
DIST.ex_amp  = 0.030; % N at severity 1.0
DIST.ex_rand = 0.015; % N rms stochastic component

% --- Gust parameters (Mars only) ---
DIST.g_rms   = 0.040; % N rms at severity 1.0
DIST.g_fc    = 1.5;   % Hz, low-pass corner of the gust spectrum

%% ---------- 0. prerequisites ----------
req = {'P','eval_ref','eval_vel','GRAV','body_names','body_g','ctrl_list', ...
       'clabels','fz_it2','fz_t1','sc_it2','sc_t1','ke_it2','ke_t1', ...
       'ki_it2','ki_t1','pids','smcs','scens'};
missing = req(~cellfun(@(v) evalin('base', sprintf('exist(''%s'',''var'')', v)) == 1, req));
if ~isempty(missing)
    error(['cdpr_dynamic_disturbance: missing workspace variables: %s\n' ...
           'Run cdpr_lowgrav_v4_4_2.m to completion first.'], strjoin(missing, ', '));
end
if any(strcmpi(DIST.modes,'gust')) && ~DIST.allow_airless
    fprintf(['[guard] ''gust'' is only physically defensible on Mars. The Moon,\n' ...
             '        Ceres and Europa have no atmosphere, so gust is evaluated\n' ...
             '        on MARS ONLY. Set DIST.allow_airless = true to override,\n' ...
             '        but do not report airless-body gust results.\n']);
end

%% ---------- 1. extract v4.4.2 functions ----------
fprintf('extracting v4.4.2 functions -> %s\n', DIST.libdir);
src = regexprep(fileread(DIST.src), '\r\n', '\n');
hdr = regexp(src, '(?m)^function [^\n]*$', 'start');
if isempty(hdr), error('no local functions found in %s', DIST.src); end
hdr(end+1) = numel(src) + 1;
if ~exist(DIST.libdir,'dir'), mkdir(DIST.libdir); end
run_sim_src = '';  n_fn = 0;
for i = 1:numel(hdr)-1
    body = src(hdr(i):hdr(i+1)-1);
    nm = regexp(body, '^function\s+(?:\[[^\]]*\]|[\w~]+)\s*=\s*(\w+)\s*\(|^function\s+(\w+)\s*\(', 'tokens','once');
    nm = nm(~cellfun(@isempty, nm));
    if isempty(nm), continue; end
    fid = fopen(fullfile(DIST.libdir,[nm{1} '.m']),'w'); fprintf(fid,'%s',body); fclose(fid);
    n_fn = n_fn + 1;
    if strcmp(nm{1},'run_sim'), run_sim_src = body; end
end
if isempty(run_sim_src), error('run_sim not found in %s', DIST.src); end
fprintf('  wrote %d functions\n', n_fn);

%% ---------- 2. build the disturbance-aware simulator ----------
% Built with strjoin, NOT sprintf: the injected block contains MATLAB
% comments and sprintf reads '%' as a format specifier, which silently
% truncates the string in MATLAB and errors in Octave.
s = strrep(run_sim_src, 'function lg = run_sim(', 'function lg = run_sim_dyn(');
NL = char(10);

% (a) state for the stochastic gust filter, initialised once
anchor_init = '    lg.qp_fail = 0;  lg.sat_steps = 0;';
assert(count_occ(s, anchor_init) == 1, 'init anchor not unique -- source changed');
init_lines = { anchor_init, ...
 '    % ---- dynamic disturbance state (study add-on) ----', ...
 '    dyn = [];  gust_state = [0;0];', ...
 '    if isfield(scen,''dyn'') && ~isempty(scen.dyn), dyn = scen.dyn; end'};
s = strrep(s, anchor_init, strjoin(init_lines, NL));

% (b) the force itself, added to the external wrench at PLANT rate so the
%     controller never sees it directly -- it is unmodelled by construction
% The disturbance is added to F_dist, which enters the plant AFTER the
% control update and is never visible to any controller -- unmodelled by
% construction, which is the whole point of the test.
anchor_w = '        Jm = jac(pos, P.anchors);';
assert(count_occ(s, anchor_w) == 1, ...
    'plant-integration anchor not found -- v4.4.2 source changed');
dyn_lines = { ...
 '        % ---- time-varying disturbance, evaluation only ----', ...
 '        if ~isempty(dyn) && dyn.sev > 0', ...
 '            switch dyn.mode', ...
 '            case ''stribeck''', ...
 '                sp = norm(vel);', ...
 '                if sp > 1e-9', ...
 '                    Fmag = dyn.Fc + (dyn.Fs - dyn.Fc)*exp(-(sp/dyn.vs)^2) ...', ...
 '                           + dyn.sig*sp;', ...
 '                    F_dist = F_dist - dyn.sev*Fmag*(vel/sp);', ...
 '                end', ...
 '            case ''extrusion''', ...
 '                F_dist = F_dist + dyn.sev*[dyn.amp*sin(2*pi*dyn.f*t); ...', ...
 '                                            dyn.amp*cos(2*pi*dyn.f*t)] ...', ...
 '                         + dyn.sev*dyn.rnd*randn(2,1);', ...
 '            case ''gust''', ...
 '                a = 1 - exp(-2*pi*dyn.fc*P.dt_plant);', ...
 '                gust_state = gust_state + a*(dyn.rms*randn(2,1)/sqrt(a/2) - gust_state);', ...
 '                F_dist = F_dist + dyn.sev*gust_state;', ...
 '            end', ...
 '        end', ...
 anchor_w};
s = strrep(s, anchor_w, strjoin(dyn_lines, NL));
assert(count_occ(s, 'dyn.sev') > 0, 'disturbance block did not survive injection');
fid = fopen(fullfile(DIST.libdir,'run_sim_dyn.m'),'w'); fprintf(fid,'%s',s); fclose(fid);
addpath(DIST.libdir); rehash;
fprintf('  built run_sim_dyn (disturbance enters the plant, not the controller)\n');

%% ---------- 3-5. sweep every mode, report, export ----------
fid = fopen('v4_dyn_disturbance.csv','w');
fprintf(fid,'Mode,Body,Controller,Severity,RMSE_mm,DegradationPct\n');
SUM = struct();

for mi = 1:numel(DIST.modes)
    mode = DIST.modes{mi};
    fprintf('\n================ mode: %s ================\n', mode);

    bodies = 1:numel(body_names);
    if strcmpi(mode,'gust') && ~DIST.allow_airless
        bodies = find(strcmpi(body_names,'Mars'));
        if isempty(bodies)
            fprintf('  Mars not in body list -- skipping gust.\n'); continue;
        end
    end
    nL = numel(DIST.levels);  nC = numel(ctrl_list);  nB = numel(bodies);
    stoch = ~strcmpi(mode,'stribeck');     % friction is deterministic
    nS = 1; if stoch, nS = DIST.n_seeds; end
    RM = nan(nB, nC, nL, nS);

    fprintf('  %d bodies x %d controllers x %d levels x %d seed(s)\n', nB, nC, nL, nS);
    t0 = tic;
    for bi = 1:nB
        b = bodies(bi);
        fprintf('  %-7s ', body_names{b});
        for c = 1:nC
            switch ctrl_list{c}
            case 'IT2', asset = struct('kind','fuzzy','fz',fz_it2,'scale',sc_it2(b),'Ke',ke_it2(b),'Ki',ki_it2(b));
            case 'T1',  asset = struct('kind','fuzzy','fz',fz_t1, 'scale',sc_t1(b), 'Ke',ke_t1(b), 'Ki',ki_t1(b));
            case 'PID', asset = pids{b};
            case 'SMC', asset = smcs{b};
            end
            for L = 1:nL
                for sd = 1:nS
                    % common random numbers: identical noise realisation
                    % across controllers at a given (level, seed), so the
                    % controller comparison is paired.
                    rng(4000 + 97*L + sd, 'twister');
                    scen = scens{1};                 % nominal plant, clean
                    scen.dyn = struct('mode',lower(mode),'sev',DIST.levels(L), ...
                        'Fc',DIST.Fc,'Fs',DIST.Fs,'vs',DIST.vs,'sig',DIST.sig, ...
                        'f',DIST.ex_f,'amp',DIST.ex_amp,'rnd',DIST.ex_rand, ...
                        'rms',DIST.g_rms,'fc',DIST.g_fc);
                    lg = run_sim_dyn(ctrl_list{c}, asset, P, body_g(b), ...
                                     eval_ref, eval_vel, scen);
                    em = sqrt(sum(lg.err.^2,2))*1000; em = em(isfinite(em));
                    if ~isempty(em), RM(bi,c,L,sd) = sqrt(mean(em.^2)); end
                end
            end
            fprintf('.');
        end
        fprintf(' done\n');
    end
    fprintf('  %.1f min\n', toc(t0)/60);

    M = mean(RM, 4, 'omitnan');
    fprintf('\n  RMSE (mm) vs severity\n');
    for bi = 1:nB
        fprintf('  %s\n', body_names{bodies(bi)});
        fprintf('    %-6s','sev');
        for c = 1:nC, fprintf('%10s', clabels{c}); end
        fprintf('%12s\n','best');
        for L = 1:nL
            fprintf('    %-6.2f', DIST.levels(L));
            for c = 1:nC, fprintf('%10.4f', M(bi,c,L)); end
            [~,ib] = min(M(bi,:,L)); fprintf('%12s\n', clabels{ib});
        end
    end

    deg = nan(nB,nC);
    for bi = 1:nB
        for c = 1:nC
            deg(bi,c) = (M(bi,c,end)-M(bi,c,1))/max(M(bi,c,1),eps)*100;
        end
    end
    fprintf('\n  degradation from severity 0 to %.1f (%%)\n', DIST.levels(end));
    fprintf('    %-8s','body');
    for c = 1:nC, fprintf('%10s', clabels{c}); end
    fprintf('\n');
    for bi = 1:nB
        fprintf('    %-8s', body_names{bodies(bi)});
        for c = 1:nC, fprintf('%9.1f%%', deg(bi,c)); end
        fprintf('\n');
    end

    % ---- verdict against the pre-registered prediction ----
    mdeg = mean(deg,1,'omitnan');
    [~,worst] = max(mdeg); [~,best] = min(mdeg);
    iSMC = find(strcmpi(ctrl_list,'SMC'));
    iFZ  = find(strcmpi(ctrl_list,'IT2') | strcmpi(ctrl_list,'T1'));
    fprintf('\n  PREDICTED: SMC degrades most, fuzzy least.\n');
    fprintf('  OBSERVED : worst = %s (%+.1f%%), best = %s (%+.1f%%)\n', ...
            clabels{worst}, mdeg(worst), clabels{best}, mdeg(best));
    if worst == iSMC && any(best == iFZ)
        fprintf('  -> CONFIRMED. Mechanism to report: inside the boundary layer the\n');
        fprintf('     SMC applies effective gain K/phi = %.0f N per unit surface,\n', ...
                smcs{bodies(1)}.K/smcs{bodies(1)}.phi);
        fprintf('     which amplifies a low-velocity disturbance; the fuzzy control\n');
        fprintf('     surface is bounded by Ku and cannot.\n');
    elseif worst == iSMC
        fprintf('  -> PARTIAL: SMC worst, but a fuzzy controller is not best.\n');
    else
        fprintf('  -> NOT CONFIRMED. Report this as it stands. Do NOT re-tune the\n');
        fprintf('     disturbance parameters to change the outcome.\n');
    end
    crossed = false;
    for bi = 1:nB
        for L = 1:nL
            if min(M(bi,iFZ,L)) < M(bi,iSMC,L)
                crossed = true;
                fprintf('  CROSSOVER on %s at severity %.2f: fuzzy %.4f < SMC %.4f mm\n', ...
                    body_names{bodies(bi)}, DIST.levels(L), min(M(bi,iFZ,L)), M(bi,iSMC,L));
            end
        end
    end
    if ~crossed
        fprintf('  No absolute crossover: SMC stays lowest at every severity.\n');
        fprintf('  A shallower degradation slope is still a reportable result, but it\n');
        fprintf('  is a DIFFERENT claim from "fuzzy beats SMC" -- do not conflate them.\n');
    end

    for bi = 1:nB
        for c = 1:nC
            for L = 1:nL
                fprintf(fid,'%s,%s,%s,%.2f,%.5f,%.2f\n', mode, body_names{bodies(bi)}, ...
                        clabels{c}, DIST.levels(L), M(bi,c,L), deg(bi,c));
            end
        end
    end
    SUM.(mode) = struct('M',M,'deg',deg,'bodies',bodies,'crossed',crossed);

    % ---- per-mode figure ----
    fw = mm2in(84); fh = mm2in(65);
    fh1 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
                 'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
    ax = axes(fh1); hold(ax,'on');
    styles = {'-o','--s',':^','-.d'};
    for c = 1:nC
        plot(ax, DIST.levels, squeeze(mean(M(:,c,:),1,'omitnan')), styles{min(c,4)}, ...
            'Color', CC{c}, 'LineWidth', 1.1, 'MarkerSize', 4, 'MarkerFaceColor','w');
    end
    lg = legend(ax, clabels, 'Location','northwest','Box','off'); fix_legend(lg);
    springer_ax(ax, sprintf('Unmodelled disturbance: %s', mode), ...
                'Disturbance severity', 'RMSE (mm)');
    print(fh1,'-dpng','-r600',sprintf('fig11_%s.png', mode));
    fprintf('  wrote fig11_%s.png\n', mode);
end
fclose(fid);

% ---- combined degradation summary across modes ----
mn = fieldnames(SUM);
if numel(mn) > 1
    Dg = zeros(numel(mn), numel(ctrl_list));
    for i = 1:numel(mn), Dg(i,:) = mean(SUM.(mn{i}).deg, 1, 'omitnan'); end
    fw = mm2in(120); fh = mm2in(70);
    f12 = figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
                 'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
    ax = axes(f12); hb = bar(ax, Dg, 'grouped');
    for c = 1:numel(ctrl_list), hb(c).FaceColor = CC{c}; hb(c).EdgeColor = 'none'; end
    set(ax,'XTick',1:numel(mn),'XTickLabel',mn);
    lg = legend(ax, clabels, 'Location','northoutside','NumColumns',2,'Box','off');
    fix_legend(lg);
    springer_ax(ax, 'Degradation Under Unmodelled Disturbance', ...
                'Disturbance class', '\DeltaRMSE (%)');
    print(f12,'-dpng','-r600','fig12_dyn_summary.png');
    fprintf('\nwrote fig12_dyn_summary.png\n');
end
fprintf('wrote v4_dyn_disturbance.csv (all modes)\n');

%% ------------------------------------------------------------------
function n = count_occ(hay, needle)
    n = numel(strfind(hay, needle));
end
