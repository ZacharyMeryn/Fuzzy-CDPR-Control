%% ========================================================================
%  run_all.m  --  execute the full study in the required order
%
%  The two follow-on studies are NOT standalone. Each reuses the workspace
%  left by the main script and reads its source file to extract shared
%  functions, so the main script must run first, to completion, in the same
%  session. This file enforces that and fails early with a useful message
%  rather than part-way through with a missing-variable error.
%
%  Usage:  cd code; run_all
%% ========================================================================

fprintf('\n=========================================================\n');
fprintf(' CDPR variable-gravity control study -- full run\n');
fprintf('=========================================================\n\n');

%% ---------- preflight ----------
fail = false;

% 1. quadprog. This is the one hard toolbox dependency and the most common
%    reason a clone of this repository does not run.
if exist('quadprog','file') ~= 2
    fprintf(2, '[FAIL] quadprog not found. The tension allocation requires the\n');
    fprintf(2, '       Optimization Toolbox. Everything in results/ is plain CSV\n');
    fprintf(2, '       and PNG if you only need to inspect the published outputs.\n');
    fail = true;
else
    fprintf('[ok]   quadprog available\n');
end

% 2. source files present, and named as the study scripts expect. Both
%    follow-on scripts locate the main file BY NAME and parse it as text.
need = {'cdpr_lowgrav_v4_4_2.m','cdpr_dynamic_disturbance.m','cdpr_fou_noise_study.m'};
for i = 1:numel(need)
    if exist(need{i},'file') ~= 2
        fprintf(2,'[FAIL] missing %s (must be on the path, and NOT renamed)\n', need{i});
        fail = true;
    else
        fprintf('[ok]   %s\n', need{i});
    end
end

% 3. clean output directory. Outputs land in pwd; anything a run does not
%    regenerate stays behind looking current, which is how a stale sweep
%    from an earlier version once came close to being cited.
stale = [dir('v4_*.csv'); dir('v4_*.tex'); dir('fig*.png')];
if ~isempty(stale)
    fprintf('\n[warn] %d pre-existing output file(s) in %s:\n', numel(stale), pwd);
    for i = 1:numel(stale)
        fprintf('         %-34s %s\n', stale(i).name, ...
                datestr(stale(i).datenum,'yyyy-mm-dd HH:MM'));
    end
    fprintf(['[warn] files this run does not regenerate will remain and will look\n' ...
             '       current. Strongly consider an empty working directory.\n']);
    resp = input('       continue anyway? [y/N] ','s');
    if isempty(resp) || lower(resp(1)) ~= 'y'
        fprintf('aborted.\n'); return;
    end
end

if fail
    error('run_all: preflight failed, see messages above.');
end

%% ---------- stage 1: main study ----------
fprintf('\n--- stage 1/3: main study (~15 min) ---\n');
t_all = tic;
cdpr_lowgrav_v4_4_2;

% the follow-on scripts need these in the BASE workspace
guard = {'P','eval_ref','eval_vel','ctrl_list','pids','smcs','ki_it2'};
if any(cellfun(@(v) evalin('base', sprintf('exist(''%s'',''var'')', v)) ~= 1, guard))
    error(['run_all: the main study did not leave the expected workspace. ' ...
           'It must complete without error before the follow-on studies run.']);
end
fprintf('\n[ok] stage 1 complete, workspace verified\n');

%% ---------- stage 2: unmodelled dynamic disturbances ----------
fprintf('\n--- stage 2/3: unmodelled time-varying disturbances ---\n');
try
    cdpr_dynamic_disturbance;
    fprintf('\n[ok] stage 2 complete\n');
catch e
    fprintf(2,'\n[FAIL] stage 2: %s\n', e.message);
    fprintf(2,'       stage 1 outputs are unaffected and remain valid.\n');
end

%% ---------- stage 3: footprint width under measurement noise ----------
fprintf('\n--- stage 3/3: FOU width under measurement noise ---\n');
try
    cdpr_fou_noise_study;
    fprintf('\n[ok] stage 3 complete\n');
catch e
    fprintf(2,'\n[FAIL] stage 3: %s\n', e.message);
    fprintf(2,'       stages 1-2 outputs are unaffected and remain valid.\n');
end

%% ---------- manifest ----------
fprintf('\n=========================================================\n');
fprintf(' total wall time: %.1f min\n', toc(t_all)/60);
made = [dir('v4_*.csv'); dir('v4_*.tex'); dir('fig*.png')];
fprintf(' outputs in %s\n', pwd);
for i = 1:numel(made)
    age = (now - made(i).datenum)*24*60;
    if age < 180, tag = 'this run'; else, tag = 'STALE -- NOT from this run'; end
    fprintf('   %-34s %s  %s\n', made(i).name, ...
            datestr(made(i).datenum,'yyyy-mm-dd HH:MM'), tag);
end
fprintf('=========================================================\n');
