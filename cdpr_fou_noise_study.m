%% ========================================================================
%  cdpr_fou_noise_study.m  --  Does the footprint of uncertainty pay under
%  measurement noise?
%
%  MOTIVATION
%  The v4.3 FOU sweep (Fig. 9) found tracking error growing monotonically
%  with FOU half-width delta under BOTH nominal and parametrically
%  mismatched conditions, and the mismatch benefit invariant to delta.
%  That plant, however, exercises only STATIC parametric uncertainty with
%  perfect state feedback. IT2 fuzzy systems are theorized to earn their
%  cost under uncertainty the T1 membership functions cannot represent:
%  measurement noise, and time-varying / state-dependent disturbance. This
%  study adds the first of those and asks whether the ranking inverts.
%
%  The hypothesis is directional and stated in advance: if the FOU does
%  anything useful, the type-reduced control surface -- flatter near the
%  origin, which COSTS precision against a constant bias -- should also be
%  LESS reactive to noise on the error signal, so wider delta should
%  degrade more slowly as noise rises. The prediction is therefore an
%  interaction (crossing curves), not a main effect.
%
%  PREREQUISITES
%  Run cdpr_lowgrav_v4_3.m to completion first, then run this file. It
%  reuses that workspace (P, eval_ref, eval_vel, tune_ref, tune_vel, GRAV,
%  ke_it2, grid_coarse, grid_fine_step, scens) and auto-extracts the v4.3
%  local functions to a folder on the path, so nothing is duplicated here
%  and the two files cannot drift apart.
%
%  OUTPUTS
%    v4_fou_noise.csv          full cell means and paired differences
%    fig9_fou_noise.png        RMSE vs noise level, one curve per FOU width
%    fig10_fou_noise_delta.png paired advantage of delta>0 over delta=0
%
%  Zach Meryn -- companion to cdpr_lowgrav_v4_3.m
%% ========================================================================

clear NOISE_CFG;
NOISE_CFG.fous        = [0 0.15 0.30 0.45];   % delta = 0 is the T1 limit
NOISE_CFG.n_seeds     = 8;                    % paired Monte Carlo replicates
NOISE_CFG.retune      = true;                 % re-tune Ku per FOU (noiseless)
NOISE_CFG.fast        = false;                % true -> 3 seeds, coarser
NOISE_CFG.src         = 'cdpr_lowgrav_v4_4_2.m';
NOISE_CFG.libdir      = fullfile(pwd, 'v43_lib');

% Noise levels. Position noise represents encoder quantization plus readout
% error; velocity noise represents the differentiation of that signal, which
% is why sigma_v is set an order of magnitude larger in relative terms --
% numerical differentiation of a quantized position at 100 Hz amplifies
% quantization by roughly 1/dt_ctrl. Levels span "good industrial encoder"
% to "degraded / dusty optical readout".
NOISE_CFG.levels = { ...
  struct('name','none',   'quant',0,     'sig_p',0,     'sig_v',0),      ...
  struct('name','low',    'quant',1e-5,  'sig_p',5e-6,  'sig_v',1.0e-3), ...
  struct('name','medium', 'quant',2e-5,  'sig_p',2e-5,  'sig_v',4.0e-3), ...
  struct('name','high',   'quant',5e-5,  'sig_p',5e-5,  'sig_v',1.0e-2), ...
  struct('name','severe', 'quant',1e-4,  'sig_p',1e-4,  'sig_v',2.5e-2)  };

if NOISE_CFG.fast, NOISE_CFG.n_seeds = 3; end

%% ---------- 0. sanity: required workspace variables ----------
req = {'P','eval_ref','eval_vel','tune_ref','tune_vel','GRAV','ke_it2', ...
       'grid_coarse','grid_fine_step'};
missing = req(~cellfun(@(v) evalin('base', sprintf('exist(''%s'',''var'')', v)) == 1, req));
if ~isempty(missing)
    error(['cdpr_fou_noise_study: missing workspace variables: %s\n' ...
           'Run cdpr_lowgrav_v4_3.m to completion first.'], strjoin(missing, ', '));
end

%% ---------- 1. extract v4.3 local functions into a callable library ----------
% Local functions inside a script are not visible to other files, so they
% are written out once as individual .m files and added to the path.
fprintf('extracting v4.3 functions -> %s\n', NOISE_CFG.libdir);
src = fileread(NOISE_CFG.src);
src = regexprep(src, '\r\n', '\n');
hdr = regexp(src, '(?m)^function [^\n]*$', 'start');
if isempty(hdr), error('no local functions found in %s', NOISE_CFG.src); end
hdr(end+1) = numel(src) + 1;
if ~exist(NOISE_CFG.libdir, 'dir'), mkdir(NOISE_CFG.libdir); end
n_fn = 0;  run_sim_src = '';
for i = 1:numel(hdr)-1
    body = src(hdr(i):hdr(i+1)-1);
    nm = regexp(body, '^function\s+(?:\[[^\]]*\]|[\w~]+)\s*=\s*(\w+)\s*\(|^function\s+(\w+)\s*\(', 'tokens', 'once');
    nm = nm(~cellfun(@isempty, nm));
    if isempty(nm), continue; end
    nm = nm{1};
    fid = fopen(fullfile(NOISE_CFG.libdir, [nm '.m']), 'w');
    fprintf(fid, '%s', body);  fclose(fid);
    n_fn = n_fn + 1;
    if strcmp(nm, 'run_sim'), run_sim_src = body; end
end
if isempty(run_sim_src), error('run_sim not found in %s', NOISE_CFG.src); end
fprintf('  wrote %d functions\n', n_fn);

%% ---------- 2. build the noise-aware run_sim variant ----------
% Four surgical substitutions on the extracted source. The TRUE error
% (e_pos, e_vel) is left untouched so that lg.err and every downstream
% metric still measure real tracking error; only the signal the CONTROLLER
% sees is corrupted. Conflating the two would make noise look like an
% improvement whenever it happened to mask error.
s = run_sim_src;
s = strrep(s, 'function lg = run_sim(', 'function lg = run_sim_noisy(');

anchor = sprintf('        if mod(k-1, P.decim) == 0\n');
assert(count_occurrences(s, anchor) == 1, 'control-update anchor not unique');
% NOTE: the injected block is built with strjoin, NOT sprintf. The block
% contains MATLAB comments, and sprintf reads '%' as a format specifier --
% it truncates the string at the first comment marker (MATLAB) or errors
% (Octave), which silently produced a run_sim_noisy with no e_pos_m in it.
NL = char(10);
inject_lines = { ...
 '            % ---- measurement noise (study add-on) ----', ...
 '            % measured = true + corruption, so the measured ERROR is the', ...
 '            % true error MINUS the corruption. Quantization applies to the', ...
 '            % position measurement itself, not to the error.', ...
 '            e_pos_m = e_pos;  e_vel_m = e_vel;', ...
 '            if isfield(scen,''noise'') && scen.noise.on', ...
 '                nq = scen.noise.quant;', ...
 '                if nq > 0', ...
 '                    e_pos_m = e_pos - (round(pos/nq)*nq - pos);', ...
 '                end', ...
 '                e_pos_m = e_pos_m - scen.noise.sig_p*randn(2,1);', ...
 '                e_vel_m = e_vel   - scen.noise.sig_v*randn(2,1);', ...
 '            end'};
inject = [anchor strjoin(inject_lines, NL) NL];
s = strrep(s, anchor, inject);

reps = { ...
 'fuzzy_ctrl(asset.fz, e_pos, e_vel, asset.scale, asset.Ke)', ...
 'fuzzy_ctrl(asset.fz, e_pos_m, e_vel_m, asset.scale, asset.Ke)', 1; ...
 'pid_step(pid, e_pos, e_vel)', 'pid_step(pid, e_pos_m, e_vel_m)', 1; ...
 's_smc  = e_vel + asset.lam*e_pos;', 's_smc  = e_vel_m + asset.lam*e_pos_m;', 1; ...
 'e_pos*P.dt_ctrl', 'e_pos_m*P.dt_ctrl', 2};
for r = 1:size(reps,1)
    n_exp = reps{r,3};
    assert(count_occurrences(s, reps{r,1}) == n_exp, ...
        'substitution %d expected %d occurrence(s) -- source changed', r, n_exp);
    s = strrep(s, reps{r,1}, reps{r,2});
end
% the generated file must actually define what it uses
assert(count_occurrences(s, 'e_pos_m = e_pos;') == 1, ...
    'noise block did not survive injection -- check for sprintf mangling');
fid = fopen(fullfile(NOISE_CFG.libdir, 'run_sim_noisy.m'), 'w');
fprintf(fid, '%s', s);  fclose(fid);
addpath(NOISE_CFG.libdir);
rehash;
fprintf('  built run_sim_noisy (true error preserved for metrics)\n');

%% ---------- 3. tune Ku per FOU width on the NOISELESS plant ----------
% Tuning is deliberately done without noise, by the identical procedure of
% Section 4.3. Letting each width tune against its own noise realization
% would confound the comparison with a tuning advantage -- the exact
% failure mode the chapter's protocol exists to prevent.
nF = numel(NOISE_CFG.fous);
FZ  = cell(1,nF);  KU = zeros(1,nF);  KI = zeros(1,nF);
fprintf('tuning Ku per FOU width (noiseless): ');
for q = 1:nF
    FZ{q} = build_fls(NOISE_CFG.fous(q));
    if NOISE_CFG.retune
        KU(q) = tune_ku(FZ{q}, P, GRAV.moon, tune_ref, tune_vel, ...
                        grid_coarse, grid_fine_step, ke_it2(1));
    else
        KU(q) = evalin('base','sc_it2(1)');
    end
    % v4.4: carry the integral gain too, tuned noiselessly by the same
    % staged procedure, so each width is evaluated in the configuration it
    % would actually be deployed in.
    KI(q) = 0;
    if evalin('base','exist(''CFG'',''var'')') && evalin('base','CFG.integral_aug')
        base_q = struct('kind','fuzzy','fz',FZ{q},'scale',KU(q),'Ke',ke_it2(1),'Ki',0);
        KI(q)  = tune_ki('IT2', base_q, P, GRAV.moon, tune_ref, tune_vel, ...
                         evalin('base','CFG.ki_grid'));
    end
    fprintf('%.2f ', KU(q));
end
fprintf('\n');

%% ---------- 4. paired Monte Carlo sweep ----------
% Common random numbers: seed is reset identically before each FOU width at
% a given (noise level, replicate), so every width sees the SAME noise
% realization. This is a paired design -- the width-to-width difference is
% computed within a shared realization, removing between-realization
% variance and making small effects detectable with few replicates.
nL = numel(NOISE_CFG.levels);
nS = NOISE_CFG.n_seeds;
RMSE  = nan(nF, nL, nS);      % full-trajectory RMSE (mm)
RHOLD = nan(nF, nL, nS);      % station-keeping RMSE (mm)
EFFRT = nan(nF, nL, nS);      % command activity, N/s (noise -> actuator work)

base_scen = struct('name','noise_eval','inject',false,'mass_f',1.00,'stiff_f',1.00);
fprintf('sweep: %d FOU x %d noise levels x %d seeds = %d runs\n', ...
        nF, nL, nS, nF*nL*nS);
t0 = tic;
for L = 1:nL
    lv = NOISE_CFG.levels{L};
    fprintf('  %-7s ', lv.name);
    for sd = 1:nS
        for q = 1:nF
            rng(1000*L + sd, 'twister');          % <-- identical across q
            scen = base_scen;
            scen.noise = struct('on', ~strcmp(lv.name,'none'), ...
                                'quant',lv.quant,'sig_p',lv.sig_p,'sig_v',lv.sig_v);
            aq = struct('kind','fuzzy','fz',FZ{q},'scale',KU(q),'Ke',ke_it2(1),'Ki',KI(q));
            lg = run_sim_noisy('IT2', aq, P, GRAV.moon, eval_ref, eval_vel, scen);

            em = sqrt(sum(lg.err.^2,2))*1000;  em = em(isfinite(em));
            if isempty(em), continue; end
            RMSE(q,L,sd) = sqrt(mean(em.^2));

            dref = sqrt(sum(diff(eval_ref).^2,2));
            kmv  = find(dref > 1e-12, 1, 'last');  if isempty(kmv), kmv = 1; end
            kh   = min(numel(em), kmv + round(0.5/P.dt_plant));
            RHOLD(q,L,sd) = sqrt(mean(em(kh:end).^2));

            Fc = lg.F(1:P.decim:end, :);
            dF = diff(Fc,1,1);
            EFFRT(q,L,sd) = mean(sqrt(sum(dF.^2,2)))/P.dt_ctrl;
        end
        fprintf('.');
    end
    fprintf(' done\n');
end
fprintf('sweep complete in %.1f min\n', toc(t0)/60);

%% ---------- 5. statistics: paired difference vs the T1 limit ----------
mR  = mean(RMSE, 3, 'omitnan');      sR = std(RMSE, 0, 3, 'omitnan');
mH  = mean(RHOLD,3, 'omitnan');
mE  = mean(EFFRT,3, 'omitnan');
DIFF = nan(nF, nL);  DSE = nan(nF, nL);  TSTAT = nan(nF, nL);
for L = 1:nL
    for q = 1:nF
        d = squeeze(RMSE(q,L,:) - RMSE(1,L,:));   % paired vs delta = 0
        d = d(isfinite(d));
        if numel(d) > 1
            DIFF(q,L)  = mean(d);
            DSE(q,L)   = std(d)/sqrt(numel(d));
            TSTAT(q,L) = mean(d)/max(std(d)/sqrt(numel(d)), eps);
        end
    end
end

fprintf('\n=== FOU x NOISE: mean RMSE (mm), %d paired seeds ===\n', nS);
fprintf('%-8s', 'noise');
for q = 1:nF, fprintf('  d=%-6.2f', NOISE_CFG.fous(q)); end
fprintf('   | best\n');
for L = 1:nL
    fprintf('%-8s', NOISE_CFG.levels{L}.name);
    for q = 1:nF, fprintf('  %7.3f ', mR(q,L)); end
    [~, ib] = min(mR(:,L));
    fprintf('   | d=%.2f\n', NOISE_CFG.fous(ib));
end

fprintf('\n=== paired advantage over the T1 limit (negative = FOU helps) ===\n');
fprintf('%-8s', 'noise');
for q = 2:nF, fprintf('   d=%-4.2f (t)   ', NOISE_CFG.fous(q)); end
fprintf('\n');
for L = 1:nL
    fprintf('%-8s', NOISE_CFG.levels{L}.name);
    for q = 2:nF
        fprintf('  %+7.4f (%+5.1f)', DIFF(q,L), TSTAT(q,L));
    end
    fprintf('\n');
end

% verdict: does any width beat delta = 0 at any noise level, with |t| > 2?
win = (DIFF(2:end,:) < 0) & (abs(TSTAT(2:end,:)) > 2);
fprintf('\nVERDICT: ');
if any(win(:))
    [qq, LL] = find(win);
    fprintf('FOU pays at %d cell(s). Strongest: ', numel(qq));
    [~, k] = min(arrayfun(@(i) DIFF(qq(i)+1, LL(i)), 1:numel(qq)));
    fprintf('delta=%.2f, noise=%s (%.4f mm, t=%.1f)\n', ...
        NOISE_CFG.fous(qq(k)+1), NOISE_CFG.levels{LL(k)}.name, ...
        DIFF(qq(k)+1,LL(k)), TSTAT(qq(k)+1,LL(k)));
    fprintf('  -> positive result: report the interaction as a design guideline.\n');
else
    fprintf('no FOU width beats the T1 limit at any tested noise level (|t|>2).\n');
    fprintf('  -> the v4.3 negative result extends to measurement noise;\n');
    fprintf('     report as a broader, better-supported null.\n');
end

%% ---------- 6. export ----------
fid = fopen('v4_fou_noise.csv','w');
fprintf(fid, ['NoiseLevel,quant_m,sig_p_m,sig_v_mps,FOU,Ku,' ...
              'RMSE_mean_mm,RMSE_std_mm,RMSEhold_mean_mm,dFrate_mean_Nps,' ...
              'PairedDiffVsT1_mm,PairedSE_mm,tstat\n']);
for L = 1:nL
    lv = NOISE_CFG.levels{L};
    for q = 1:nF
        fprintf(fid, '%s,%.3e,%.3e,%.3e,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f,%.5f,%.2f\n', ...
            lv.name, lv.quant, lv.sig_p, lv.sig_v, NOISE_CFG.fous(q), KU(q), ...
            mR(q,L), sR(q,L), mH(q,L), mE(q,L), DIFF(q,L), DSE(q,L), TSTAT(q,L));
    end
end
fclose(fid);
fprintf('\nwrote v4_fou_noise.csv\n');

%% ---------- 7. figures ----------
xs = 1:nL;  lbl = cellfun(@(c) c.name, NOISE_CFG.levels, 'UniformOutput', false);
styles = {'-o','--s',':^','-.d'};

f9 = figure('Color','w','Position',[100 100 560 380]);
ax = axes(f9); hold(ax,'on');
for q = 1:nF
    errorbar(ax, xs, mR(q,:), sR(q,:)/sqrt(nS), styles{min(q,4)}, ...
        'LineWidth',1.4, 'Color', [1 1 1]*0.22*(q-1), 'MarkerFaceColor','w');
end
set(ax,'XTick',xs,'XTickLabel',lbl,'XLim',[0.7 nL+0.3]);
legend(ax, arrayfun(@(v) sprintf('\\delta = %.2f', v), NOISE_CFG.fous, ...
    'UniformOutput', false), 'Location','northwest','Box','off');
try
    springer_ax(ax, 'Effect of FOU Width Under Measurement Noise', ...
        'Noise level', 'RMSE (mm)');
catch
    title(ax,'Effect of FOU Width Under Measurement Noise');
    xlabel(ax,'Noise level'); ylabel(ax,'RMSE (mm)'); grid(ax,'on');
end
print(f9,'-dpng','-r600','fig9_fou_noise.png');

f10 = figure('Color','w','Position',[100 100 560 380]);
ax = axes(f10); hold(ax,'on');
for q = 2:nF
    errorbar(ax, xs, DIFF(q,:), DSE(q,:), styles{min(q,4)}, ...
        'LineWidth',1.4, 'Color', [1 1 1]*0.22*(q-1), 'MarkerFaceColor','w');
end
yline(ax, 0, 'k-', 'LineWidth', 1.0);
set(ax,'XTick',xs,'XTickLabel',lbl,'XLim',[0.7 nL+0.3]);
legend(ax, arrayfun(@(v) sprintf('\\delta = %.2f', v), NOISE_CFG.fous(2:end), ...
    'UniformOutput', false), 'Location','northwest','Box','off');
try
    springer_ax(ax, 'Paired Advantage Over the Type-1 Limit', ...
        'Noise level', '\DeltaRMSE vs \delta = 0 (mm)');
catch
    title(ax,'Paired Advantage Over the Type-1 Limit');
    xlabel(ax,'Noise level'); ylabel(ax,'\DeltaRMSE vs \delta = 0 (mm)'); grid(ax,'on');
end
print(f10,'-dpng','-r600','fig10_fou_noise_delta.png');
fprintf('wrote fig9_fou_noise.png, fig10_fou_noise_delta.png\n');

%% ------------------------------------------------------------------
function n = count_occurrences(hay, needle)
    n = numel(strfind(hay, needle));
end
