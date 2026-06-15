% impactful asf  —  V3 (fair four-way comparison: robust tuning for all)
% four-way controller showdown — IT2 vs T3 vs FOPID vs neuro, reduced gravity
% CDPR 1m x 1m, four planetary bodies
%
% same plant / scenarios / springer figure style as LowGrav_4_21_2pm.m.
% IT2 is NOT modified (verbatim).  goal of v3: give every controller an
% honest, well-tuned shot so the comparison is clean.
%
% ---- WHAT CHANGED FROM V2 (and why) -------------------------------------
% 1. Robust scale tuning. v2 tuned the IT2/T3 output scale by finite-difference
%    gradient descent, which was UNSTABLE on the T3 surface (its cost curve
%    spiked and never converged, so T3 "blew up" to ~220 mm on Ceres/Europa/
%    Mars — a tuning artifact, not a property of Type-3). v3 replaces it with a
%    bullet-proof 1-D grid scan + local refine: the output scale is a single
%    number, so a grid search cannot oscillate and finds the true optimum for
%    BOTH IT2 and T3. This gives T3 a fair shot.
% 2. Scale ceiling reverted 5.0 -> 3.0. The v2 raise (3->5) over-drove the
%    fuzzy controllers at higher gravity and made Mars IT2 oscillate (172 mm
%    vs the paper's 30 mm). 3.0 is the original, better-calibrated bound.
% 3. The learning-curve figure is replaced by a cost-vs-scale tuning curve
%    (Moon), which actually visualises that the tuning is well-posed.
% FOPID and the NN/ANFIS paths are unchanged from v2 (they already work).
%
% !!! NOT EXECUTED IN THE AUTHORING ENVIRONMENT — review the printout. !!!
% needs: Optimization Toolbox (quadprog), Fuzzy Logic Toolbox (mamfis/evalfis;
%        genfis/anfis for the NN — ffnn backprop fallback used otherwise)

clear; clc; close all;
warning('off','fuzzy:general:warnEvalfis_OutOfRangeInput');


% =========================================================================
% BLOCK 1: PARAMETERS  (identical to LowGrav_4_21_2pm.m)
% =========================================================================
GRAVITY.earth  = 9.81;
GRAVITY.moon   = 1.62;
GRAVITY.mars   = 3.72;
GRAVITY.ceres  = 0.27;
GRAVITY.europa = 1.32;

P.anchors = [0.0 1.0; 1.0 1.0; 1.0 0.0; 0.0 0.0];
P.num_cables = 4;
P.ws_min = [0.05 0.05];  P.ws_max = [0.95 0.95];
P.mass_EE = 0.5;
P.tube_force_per_meter = 0.3;
P.cable_diameter = 0.001;  P.cable_density = 0.97;
P.cable_mass_per_m = P.cable_density * (pi*(P.cable_diameter/2)^2) * 1e6;
P.cable_stiffness = 50000;  P.cable_damping = 5.0;
P.cable_compliance_bandwidth = 20.0;
P.T_min = 2.0;  P.T_max = 50.0;
P.dt = 0.01;  P.t_total = 10.0;  P.t_vec = 0:P.dt:P.t_total;  P.N = numel(P.t_vec);
P.disturbance_time = 2.5;  P.disturbance_magnitude = 1.5;  P.disturbance_duration = 0.3;
P.vel_max = 0.20;  P.accel_max = 0.40;
P.mass_uncertainty = 0.15;  P.stiffness_uncertainty = 0.15;

fprintf('four-controller comparison V3 loaded\n');
fprintf('bodies: moon (%.2f) | mars (%.2f) | ceres (%.2f) | europa (%.2f)\n', ...
        GRAVITY.moon, GRAVITY.mars, GRAVITY.ceres, GRAVITY.europa);
fprintf('100 Hz | %.0f s | %d steps\n\n', P.t_total, P.N);

if exist('mamfis','file') ~= 2 || exist('quadprog','file') ~= 2
    error('need Fuzzy Logic Toolbox (mamfis) and Optimization Toolbox (quadprog).');
end
has_anfis = (exist('anfis','file')==2) && (exist('genfis','file')==2);
fprintf('ANFIS available: %d (ffnn backprop fallback used if 0)\n\n', has_anfis);


% =========================================================================
% BLOCK 2: build all four controllers + tune them
% =========================================================================
body_names = {'Moon','Mars','Ceres','Europa'};
body_g     = [GRAVITY.moon, GRAVITY.mars, GRAVITY.ceres, GRAVITY.europa];
nb = numel(body_names);

fprintf('building IT2 fuzzy controller...\n');
it2_fis = build_it2_fls(0.15);
fprintf('  %d rules\n', numel(it2_fis.Rules));

fprintf('building T3 fuzzy controller...\n');
t3 = build_t3_fls(0.15, 0.30, 5);
fprintf('  %d rules | %d alpha-planes\n', t3.n_rules, numel(t3.alpha));

t_ref = linspace(0,3,300);
[ref_tune,~] = make_trap_traj([0.3 0.3],[0.7 0.7], t_ref, P.vel_max, P.accel_max);
it2_scale = zeros(1,nb);  t3_scale = zeros(1,nb);
scan_it2 = cell(1,nb);  scan_t3 = cell(1,nb);
fprintf('\noutput-scale tuning (robust 1-D grid search, IT2 + T3):\n');
for b = 1:nb
    fprintf('  %-7s ', body_names{b});
    [it2_scale(b), scan_it2{b}] = tune_scale(@it2_eval_wrap, it2_fis, ref_tune, P, body_g(b));
    [t3_scale(b),  scan_t3{b}]  = tune_scale(@t3_eval_wrap,  t3,      ref_tune, P, body_g(b));
    fprintf(' IT2=%.3f | T3=%.3f\n', it2_scale(b), t3_scale(b));
end

fprintf('\nFOPID tuning (ziegler-nichols seed + fminsearch, bounded gains):\n');
[step_ref, step_vel] = make_trap_traj([0.3 0.4],[0.7 0.6], P.t_vec, P.vel_max, P.accel_max);
[Ku,Tu] = estimate_ultimate_gain(P, GRAVITY.moon, step_ref, step_vel);
if isfinite(Ku) && isfinite(Tu) && Ku>0 && Tu>0
    x0 = zn_baseline(Ku, Tu);
    fprintf('  ZN seed: Ku=%.2f Tu=%.3fs -> Kp=%.2f Ki=%.2f Kd=%.2f, lambda=mu=1\n', ...
            Ku, Tu, x0(1), x0(2), x0(3));
else
    x0 = [50 2.0 8.0 1.0 1.0];
    fprintf('  ZN inconclusive; seeding from bandwidth PID [50 2 8], lambda=mu=1\n');
end
fopid_base = fopid_tune(@(x) fopid_step_cost(x, P, GRAVITY.moon, step_ref, step_vel), x0);
fprintf('  tuned (Moon): Kp=%.3f Ki=%.3f Kd=%.3f lambda=%.3f mu=%.3f\n', ...
        fopid_base(1),fopid_base(2),fopid_base(3),fopid_base(4),fopid_base(5));

fprintf('\nNN training (teacher = IT2 reference, full-surface + trajectory data)...\n');
[Xtr, Ytr] = build_nn_dataset(it2_fis, it2_scale, body_g, P, step_ref, step_vel);
nn_model = nn_train(Xtr, Ytr);
fprintf('  trained (%s) on %d samples\n', nn_model.type, size(Xtr,1));

A.it2_fis=it2_fis; A.it2_scale=it2_scale;
A.t3=t3; A.t3_scale=t3_scale;
A.fopid_base=fopid_base;
A.nn_model=nn_model;

controllers = {'IT2','T3','FOPID','NN'};
clabels     = {'IT2-FLS','T3-FLS','FOPID','Neuro (ANFIS/FFNN)'};
nc = numel(controllers);


% =========================================================================
% BLOCK 3: RUN EVERYTHING — 4 bodies x 4 controllers x 2 tests
% =========================================================================
fprintf('\nrunning simulation...\n');
nominal = struct('mass_factor',1.0,'stiffness_factor',1.0);
results = struct();

for b = 1:nb
    fprintf('\n  %s (g=%.2f)\n', body_names{b}, body_g(b));
    for k = 1:nc
        ck = controllers{k};
        t0 = tic;
        log1 = run_sim(ck, b, A, P, body_g(b), step_ref, step_vel, false, nominal);
        log2 = run_sim(ck, b, A, P, body_g(b), step_ref, step_vel, true,  nominal);
        rt = toc(t0);
        m = compute_metrics_ext(log1, log2, P);  m.runtime = rt;

        fn = sprintf('%s_%s', body_names{b}, ck);
        results.(fn).step = log1;  results.(fn).dist = log2;  results.(fn).m = m;

        stab=''; if ~m.stable, stab='  [UNSTABLE/NON-CONVERGED]'; end
        fprintf('    %-6s RMSE=%7.2fmm settle=%5.2fs ISE=%9.1f IAE=%7.1f rt=%.2fs%s\n', ...
                ck, m.rmse, m.settling, m.ISE, m.IAE, rt, stab);
    end
end


% =========================================================================
% BLOCK 4: CONSOLE TABLE
% =========================================================================
fprintf('\n');
fprintf('%-18s  %-8s  %-9s  %-8s  %-9s  %-9s  %-7s\n', ...
        'body + controller','RMSE','MaxErr','Settle','ISE','IAE','rt(s)');
fprintf('%s\n', repmat('-',1,82));
for b = 1:nb
    for k = 1:nc
        fn = sprintf('%s_%s', body_names{b}, controllers{k});  m = results.(fn).m;
        fprintf('%-18s  %6.2fmm  %7.2fmm  %6.2fs  %9.1f  %8.1f  %5.2f\n', ...
                [body_names{b} ' ' controllers{k}], m.rmse, m.max_err, m.settling, m.ISE, m.IAE, m.runtime);
    end
end


% =========================================================================
% BLOCK 5: EXPORT CSV
% =========================================================================
fid = fopen('all_controllers_results_v3.csv','w');
fprintf(fid,['Body,Controller,Stable,RMSE_mm,MaxErr_mm,Settle_s,Recovery_s,', ...
             'Overshoot_pct,SteadyStateErr_mm,ISE_mm2s,IAE_mms,Runtime_s,TViol\n']);
for b = 1:nb
    for k = 1:nc
        fn = sprintf('%s_%s', body_names{b}, controllers{k});  m = results.(fn).m;
        rec = m.recovery; if isnan(rec), recs='NaN'; else, recs=sprintf('%.3f',rec); end
        fprintf(fid,'%s,%s,%d,%.4f,%.4f,%.4f,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n', ...
            body_names{b}, clabels{k}, m.stable, m.rmse, m.max_err, m.settling, recs, ...
            m.overshoot, m.sse, m.ISE, m.IAE, m.runtime, m.tension_viol);
    end
end
fclose(fid);
fprintf('\nexported all_controllers_results_v3.csv\n');


% =========================================================================
% BLOCK 6: FIGURES — springer journal style (series = the four controllers)
% =========================================================================
mm2in = @(x) x/25.4;
FN='Times New Roman'; FS=8; LW=1.2; LWt=0.6;
CC = {[0 0 0],[0.40 0.40 0.40],[0.65 0.65 0.65],[0.15 0.15 0.15]};
LSc = {'-','--',':','-.'};
set(0,'DefaultAxesFontName',FN,'DefaultAxesFontSize',FS, ...
      'DefaultTextFontName',FN,'DefaultTextFontSize',FS, ...
      'DefaultLegendFontSize',7,'DefaultLegendFontName',FN);

RMSE=zeros(nb,nc); SETT=zeros(nb,nc); ISEm=zeros(nb,nc); IAEm=zeros(nb,nc);
for b=1:nb, for k=1:nc
    m=results.(sprintf('%s_%s',body_names{b},controllers{k})).m;
    RMSE(b,k)=m.rmse; SETT(b,k)=m.settling; ISEm(b,k)=m.ISE; IAEm(b,k)=m.IAE;
end, end

ymax=0;
for b=1:nb, for k=1:nc
    e=sqrt(sum(results.(sprintf('%s_%s',body_names{b},controllers{k})).step.error.^2,2))*1000;
    ymax=max(ymax,max(e(isfinite(e))));
end, end
y_shared = ymax*1.1;

% fig 1: RMSE + settling grouped bars
fw=mm2in(174); fh=mm2in(70);
fig1=figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax1=subplot(1,2,1); hb=bar(ax1,RMSE,'grouped'); for k=1:nc, hb(k).FaceColor=CC{k}; hb(k).EdgeColor='none'; end
set(ax1,'XTick',1:nb,'XTickLabel',body_names);
lg=legend(ax1,clabels,'Location','northwest','Box','off'); fix_legend(lg);
springer_ax(ax1,'(a) Step Response RMSE','Planetary body','RMSE (mm)');
ax2=subplot(1,2,2); hb=bar(ax2,SETT,'grouped'); for k=1:nc, hb(k).FaceColor=CC{k}; hb(k).EdgeColor='none'; end
set(ax2,'XTick',1:nb,'XTickLabel',body_names);
lg=legend(ax2,clabels,'Location','northeast','Box','off'); fix_legend(lg);
springer_ax(ax2,'(b) Settling Time','Planetary body','Time (s)');
print(fig1,'-dpng','-r600','fig1_cmp_rmse_settling.png');

% fig 2: tracking error time series, one panel per controller, shared y
fw=mm2in(174); fh=mm2in(140);
fig2=figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ptag = {'(a) ','(b) ','(c) ','(d) '};
for k=1:nc
    axk=subplot(2,2,k);
    for b=1:nb
        e=sqrt(sum(results.(sprintf('%s_%s',body_names{b},controllers{k})).step.error.^2,2))*1000;
        plot(axk,P.t_vec,e,LSc{b},'Color',CC{b},'LineWidth',LW); hold(axk,'on');
    end
    ylim(axk,[0 y_shared]);
    if k==1, lg=legend(axk,body_names,'Location','northeast','Box','off'); fix_legend(lg); end
    springer_ax(axk,[ptag{k} clabels{k}],'Time (s)','Position error (mm)');
end
print(fig2,'-dpng','-r600','fig2_cmp_tracking_error.png');

% fig 3: disturbance rejection on the Moon, all four overlaid
fw=mm2in(84); fh=mm2in(72);
fig3=figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax5=axes(fig3);
for k=1:nc
    e=sqrt(sum(results.(sprintf('Moon_%s',controllers{k})).dist.error.^2,2))*1000;
    plot(ax5,P.t_vec,e,LSc{k},'Color',CC{k},'LineWidth',LW); hold(ax5,'on');
end
xline(ax5,P.disturbance_time,':','Color',[0.5 0.5 0.5],'LineWidth',LWt,'HandleVisibility','off');
xline(ax5,P.disturbance_time+P.disturbance_duration,':','Color',[0.5 0.5 0.5],'LineWidth',LWt,'HandleVisibility','off');
lg=legend(ax5,clabels,'Location','northwest','Box','off'); fix_legend(lg);
springer_ax(ax5,'Disturbance Rejection (Moon)','Time (s)','Position error (mm)');
print(fig3,'-dpng','-r600','fig3_cmp_disturbance.png');

% fig 4: ISE + IAE grouped bars
fw=mm2in(174); fh=mm2in(70);
fig4=figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
            'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
ax6=subplot(1,2,1); hb=bar(ax6,ISEm,'grouped'); for k=1:nc, hb(k).FaceColor=CC{k}; hb(k).EdgeColor='none'; end
set(ax6,'XTick',1:nb,'XTickLabel',body_names);
lg=legend(ax6,clabels,'Location','northeast','Box','off'); fix_legend(lg);
springer_ax(ax6,'(a) Integral Squared Error','Planetary body','ISE (mm^2 s)');
ax7=subplot(1,2,2); hb=bar(ax7,IAEm,'grouped'); for k=1:nc, hb(k).FaceColor=CC{k}; hb(k).EdgeColor='none'; end
set(ax7,'XTick',1:nb,'XTickLabel',body_names);
lg=legend(ax7,clabels,'Location','northeast','Box','off'); fix_legend(lg);
springer_ax(ax7,'(b) Integral Absolute Error','Planetary body','IAE (mm s)');
print(fig4,'-dpng','-r600','fig4_cmp_ise_iae.png');

% fig 5: output-scale tuning curve (Moon) — cost vs scale for IT2 and T3
if ~isempty(scan_it2{1})
    fw=mm2in(84); fh=mm2in(65);
    fig5=figure('Color','white','Units','inches','Position',[1 1 fw fh], ...
                'PaperUnits','inches','PaperSize',[fw fh],'PaperPosition',[0 0 fw fh]);
    ax8=axes(fig5);
    s1=scan_it2{1}; s3=scan_t3{1};
    plot(ax8,s1.scales,s1.costs,'-','Color',CC{1},'LineWidth',LW); hold(ax8,'on');
    plot(ax8,s3.scales,s3.costs,'--','Color',CC{2},'LineWidth',LW);
    [~,i1]=min(s1.costs); [~,i3]=min(s3.costs);
    plot(ax8,s1.scales(i1),s1.costs(i1),'o','Color',CC{1},'MarkerSize',4,'MarkerFaceColor',CC{1},'HandleVisibility','off');
    plot(ax8,s3.scales(i3),s3.costs(i3),'s','Color',CC{2},'MarkerSize',4,'MarkerFaceColor',CC{2},'HandleVisibility','off');
    ylim(ax8,[0, 6*min([s1.costs(:);s3.costs(:)])]);
    lg=legend(ax8,{'IT2-FLS','T3-FLS'},'Location','northwest','Box','off'); fix_legend(lg);
    springer_ax(ax8,'Output-Scale Tuning (Moon)','Output scale','Tracking cost \itJ');
    print(fig5,'-dpng','-r600','fig5_cmp_scale_tuning.png');
end

fprintf('\nfigures saved at 600 DPI. done.\n');


% #########################################################################
% LOCAL FUNCTIONS
% #########################################################################

function [pos_traj, vel_traj] = make_trap_traj(pos_start, pos_end, t_vec, vel_max, accel_max)
    N=numel(t_vec); d=norm(pos_end-pos_start); dir=(pos_end-pos_start)/max(d,1e-9);
    v_peak=min(vel_max,sqrt(accel_max*d)); t_ramp=v_peak/accel_max;
    t_flat=(d-accel_max*t_ramp^2)/v_peak;
    if t_flat<0, t_ramp=sqrt(d/accel_max); v_peak=accel_max*t_ramp; t_flat=0; end
    te1=t_ramp; te2=t_ramp+t_flat; te3=2*t_ramp+t_flat;
    pos_traj=zeros(N,2); vel_traj=zeros(N,2);
    for k=1:N
        t=t_vec(k);
        if t<=te1,     s=0.5*accel_max*t^2;                          v=accel_max*t;
        elseif t<=te2, s=0.5*accel_max*t_ramp^2+v_peak*(t-te1);       v=v_peak;
        elseif t<=te3, d2=t-te2; s=0.5*accel_max*t_ramp^2+v_peak*t_flat+v_peak*d2-0.5*accel_max*d2^2;
                       v=v_peak-accel_max*d2;
        else,          s=d; v=0;
        end
        s=min(s,d); pos_traj(k,:)=pos_start+s*dir; vel_traj(k,:)=v*dir;
    end
end

function [L_arc, sag_mid, L_chord] = ik_with_sag(pos, anchors, T_est, g, cmpm)
    L_chord=sqrt(sum((anchors-pos).^2,2)); w=cmpm*g; T_safe=max(T_est,0.1);
    sag_factor=(w.*L_chord).^2 ./ (24.*T_safe.^2);
    L_arc=L_chord.*(1+sag_factor); sag_mid=w.*L_chord.^2 ./ (8.*T_safe);
end

function J = compute_jacobian(pos, anchors)
    L=sqrt(sum((anchors-pos).^2,2)); J=(anchors-pos)./L;
end

function [T, feasible] = tension_qp(pos, anchors, F_ext, T_min, T_max)
    n=size(anchors,1);
    if ~all(isfinite(F_ext(:)))
        T=(T_min+T_max)/2*ones(n,1); feasible=false; return;
    end
    J=compute_jacobian(pos,anchors);
    opts=optimoptions('quadprog','Display','off');
    [T,~,flag]=quadprog(eye(n),zeros(n,1),[],[],J',F_ext(:),T_min*ones(n,1),T_max*ones(n,1),[],opts);
    feasible=(flag==1); if ~feasible, T=(T_min+T_max)/2*ones(n,1); end
end

% ----- IT2 FLS build (VERBATIM from LowGrav_4_21_2pm.m, do not edit) -----
function fis = build_it2_fls(sigma_uncertainty)
    fis=mamfis('Name','IT2_PositionController','AndMethod','min','OrMethod','max', ...
               'ImplicationMethod','min','AggregationMethod','max','DefuzzificationMethod','centroid');
    fis=addInput(fis,[-0.4 0.4],'Name','pos_error'); sig_e=0.06; su=sig_e*(1+sigma_uncertainty);
    ce=[-0.25 -0.10 0.0 0.10 0.25]; mf={'NB','NS','ZE','PS','PB'};
    for i=1:5, fis=addMF(fis,'pos_error','gaussmf',[su ce(i)],'Name',mf{i}); end
    fis=addInput(fis,[-0.3 0.3],'Name','vel_error'); svu=0.05*(1+sigma_uncertainty);
    cv=[-0.20 -0.08 0.0 0.08 0.20];
    for i=1:5, fis=addMF(fis,'vel_error','gaussmf',[svu cv(i)],'Name',mf{i}); end
    fis=addOutput(fis,[-20 20],'Name','force_cmd'); oc=[-16 -8 0 8 16];
    for i=1:5, fis=addMF(fis,'force_cmd','trimf',[oc(i)-4 oc(i) oc(i)+4],'Name',mf{i}); end
    rt=[1 1 1 2 2;1 1 2 2 3;1 2 3 4 5;3 4 4 5 5;4 4 5 5 5]; rules=[];
    for r=1:5, for c=1:5, rules=[rules; r c rt(r,c) 1 1]; end, end %#ok<AGROW>
    fis=addRule(fis,rules);
end

function F_cmd = it2_evaluate(fis, pos_error, vel_error, output_scale)
    if isempty(fis), F_cmd=[0;0]; return; end
    ex=max(-0.4,min(0.4,pos_error(1))); ey=max(-0.4,min(0.4,pos_error(2)));
    vx=max(-0.3,min(0.3,vel_error(1))); vy=max(-0.3,min(0.3,vel_error(2)));
    Fx=evalfis(fis,[ex vx])*output_scale; Fy=evalfis(fis,[ey vy])*output_scale;
    F_cmd=[Fx;Fy];
end

% ===================== TYPE-3 FLS =====================
function t3 = build_t3_fls(sigma_fou, beta_t3, n_alpha)
    t3.e_sigma=0.06; t3.e_centers=[-0.25 -0.10 0.0 0.10 0.25];
    t3.v_sigma=0.05; t3.v_centers=[-0.20 -0.08 0.0 0.08 0.20];
    t3.out_centers=[-16 -8 0 8 16];
    t3.rule_table=[1 1 1 2 2;1 1 2 2 3;1 2 3 4 5;3 4 4 5 5;4 4 5 5 5];
    t3.sigma_fou=sigma_fou; t3.beta_t3=beta_t3;
    t3.alpha=linspace(1/n_alpha,1.0,n_alpha);
    rc=[]; ii=[]; jj=[];
    for i=1:5, for j=1:5
        o=t3.rule_table(i,j); rc=[rc;t3.out_centers(o)]; ii=[ii;i]; jj=[jj;j]; %#ok<AGROW>
    end, end
    t3.rule_centroids=rc; t3.rule_i=ii; t3.rule_j=jj; t3.n_rules=numel(rc);
end

function F_cmd = t3_evaluate(t3, pos_error, vel_error, output_scale)
    if isempty(t3), F_cmd=[0;0]; return; end
    ex=max(-0.4,min(0.4,pos_error(1))); ey=max(-0.4,min(0.4,pos_error(2)));
    vx=max(-0.3,min(0.3,vel_error(1))); vy=max(-0.3,min(0.3,vel_error(2)));
    Fx=t3_infer_axis(t3,ex,vx)*output_scale; Fy=t3_infer_axis(t3,ey,vy)*output_scale;
    F_cmd=[Fx;Fy];
end

function y = t3_infer_axis(t3, e, v)
    num=0; den=0;
    for k=1:numel(t3.alpha)
        a=t3.alpha(k);
        fou_k=t3.sigma_fou*(1 - t3.beta_t3*(1-a));
        e_su=t3.e_sigma*(1+fou_k); e_sl=t3.e_sigma*(1-fou_k);
        v_su=t3.v_sigma*(1+fou_k); v_sl=t3.v_sigma*(1-fou_k);
        mu_e_up=gaussv(e,t3.e_centers,e_su); mu_e_lo=gaussv(e,t3.e_centers,e_sl);
        mu_v_up=gaussv(v,t3.v_centers,v_su); mu_v_lo=gaussv(v,t3.v_centers,v_sl);
        w_hi=min(mu_e_up(t3.rule_i),mu_v_up(t3.rule_j));
        w_lo=min(mu_e_lo(t3.rule_i),mu_v_lo(t3.rule_j));
        [yl,yr]=km_reduce(t3.rule_centroids,w_lo,w_hi);
        num=num+a*0.5*(yl+yr); den=den+a;
    end
    if den<eps, y=0; else, y=num/den; end
end

function m = gaussv(x, centers, sigma)
    m=exp(-((x-centers).^2)./(2*sigma.^2));
end

function [yl,yr] = km_reduce(c, w_lo, w_hi)
    c=c(:); w_lo=w_lo(:); w_hi=w_hi(:);
    if sum(w_hi)<eps, yl=0; yr=0; return; end
    [c,idx]=sort(c); w_lo=w_lo(idx); w_hi=w_hi(idx); N=numel(c);
    w=0.5*(w_lo+w_hi); yl=wsum(c,w);
    for it=1:N
        L=find(c<=yl,1,'last'); if isempty(L), L=1; end, if L>=N, L=N-1; end
        wn=w_hi; wn(L+1:end)=w_lo(L+1:end); yln=wsum(c,wn);
        if abs(yln-yl)<1e-9, yl=yln; break; end, yl=yln;
    end
    w=0.5*(w_lo+w_hi); yr=wsum(c,w);
    for it=1:N
        R=find(c<=yr,1,'last'); if isempty(R), R=1; end, if R>=N, R=N-1; end
        wn=w_lo; wn(R+1:end)=w_hi(R+1:end); yrn=wsum(c,wn);
        if abs(yrn-yr)<1e-9, yr=yrn; break; end, yr=yrn;
    end
end

function y = wsum(c,w)
    s=sum(w); if s<eps, y=0; else, y=sum(c.*w)/s; end
end

% ===================== FRACTIONAL-ORDER PID =====================
function w = gl_weights(ord, Lmem)
    w=zeros(1,Lmem+1); w(1)=1;
    for j=1:Lmem, w(j+1)=(1-(ord+1)/j)*w(j); end
end

function st = fopid_init(params, dt, output_lim, Lmem)
    st.Kp=params(1); st.Ki=params(2); st.Kd=params(3);
    st.lambda=params(4); st.mu=params(5); st.dt=dt;
    st.output_lim=output_lim; st.Lmem=Lmem;
    st.w_int=gl_weights(-st.lambda,Lmem); st.w_der=gl_weights(st.mu,Lmem);
    st.e_hist=zeros(2,Lmem+1); st.n_seen=0;
end

function [F_cmd, st] = fopid_update(st, pos_des, pos_est)
    e=(pos_des(:)-pos_est(:));
    st.e_hist=[e, st.e_hist(:,1:end-1)]; st.n_seen=min(st.n_seen+1,st.Lmem+1);
    nu=st.n_seen; Eh=st.e_hist(:,1:nu);
    int_term=(st.dt^st.lambda)*(Eh*st.w_int(1:nu).');
    der_term=(st.dt^(-st.mu))*(Eh*st.w_der(1:nu).');
    F_cmd=st.Kp*e + st.Ki*int_term + st.Kd*der_term;
    fm=norm(F_cmd); if fm>st.output_lim, F_cmd=F_cmd*st.output_lim/fm; end
end

function x0 = zn_baseline(Ku, Tu)
    x0=[0.6*Ku, 1.2*Ku/Tu, 0.075*Ku*Tu, 1.0, 1.0];
end

function p = fopid_tune(cost_fn, x0)
    opts=optimset('Display','off','MaxIter',400,'MaxFunEvals',800,'TolFun',1e-4,'TolX',1e-4);
    xo=fminsearch(@(x) guarded_cost(x,cost_fn), x0, opts); p=clip_params(xo);
end

function J = guarded_cost(x, cost_fn)
    xp=clip_params(x); pen=1e3*sum((x-xp).^2); J=cost_fn(xp)+pen;
    if ~isfinite(J), J=1e9; end
end

function xp = clip_params(x)
    xp=x;
    xp(1)=min(150,max(0,x(1)));
    xp(2)=min(150,max(0,x(2)));
    xp(3)=min( 80,max(0,x(3)));
    xp(4)=min(1.5,max(0.05,x(4)));
    xp(5)=min(1.5,max(0.05,x(5)));
end

function [Ku,Tu] = estimate_ultimate_gain(P, g, ref, vel)
    Ku=NaN; Tu=NaN;
    for Kp=[5 10 20 40 80 160 320]
        Atmp.fopid_base=[Kp 0 0 1 1];
        log=run_sim('FOPID',1,Atmp,P,g,ref,vel,false,struct('mass_factor',1,'stiffness_factor',1));
        if ~log.stable, break; end
        e=log.error(end-round(4/P.dt):end,1); e=e-mean(e);
        zc=sum(abs(diff(sign(e)))>0); amp=std(e);
        if zc>=6 && amp>1e-3, Ku=Kp; Tu=2*(numel(e)*P.dt)/zc; return; end
    end
end

function J = fopid_step_cost(x, P, g, ref, vel)
    Atmp.fopid_base=x;
    log=run_sim('FOPID',1,Atmp,P,g,ref,vel,false,struct('mass_factor',1,'stiffness_factor',1));
    if ~log.stable, J=1e9; return; end
    e_mag=sqrt(sum(log.error.^2,2))*1000; ise=sum(e_mag.^2)*P.dt;
    ss=mean(e_mag(end-round(0.5/P.dt):end));
    os=max(0,-min(log.error*((ref(end,:)/norm(ref(end,:)))'))*1000);
    eff=mean(sum(log.F_cmd.^2,2));
    J=ise + 50*ss + 5*os + 0.5*eff;
    if ~isfinite(J), J=1e9; end
end

% ===================== NEURO (ANFIS / FFNN) =====================
function model = nn_train(X, Y)
    X=X(:,1:2); Y=Y(:);
    good=all(isfinite(X),2) & isfinite(Y); X=X(good,:); Y=Y(good);
    if (exist('anfis','file')==2) && (exist('genfis','file')==2)
        gopt=genfisOptions('GridPartition');
        gopt.NumMembershipFunctions=[3 3];
        gopt.InputMembershipFunctionType=["gaussmf" "gaussmf"];
        initFis=genfis(X,Y,gopt);
        aopt=anfisOptions('InitialFIS',initFis,'EpochNumber',40, ...
            'DisplayANFISInformation',0,'DisplayErrorValues',0, ...
            'DisplayStepSize',0,'DisplayFinalResults',0);
        model.type='anfis'; model.fis=anfis([X Y],aopt);
        model.in_range=[model.fis.Inputs(1).Range; model.fis.Inputs(2).Range];
    else
        model=train_ffnn(X,Y);
    end
end

function model = train_ffnn(X, Y)
    H=12; epochs=4000; lr=0.02;
    xmu=mean(X,1); xsig=std(X,0,1)+1e-9; ymu=mean(Y); ysig=std(Y)+1e-9;
    Xn=(X-xmu)./xsig; Yn=(Y-ymu)./ysig; M=size(Xn,1);
    rng(0); W1=0.5*randn(H,2); b1=zeros(H,1); W2=0.5*randn(1,H); b2=0;
    Xt=Xn.';
    for ep=1:epochs
        Z1=W1*Xt+b1; A1=tanh(Z1); Yh=W2*A1+b2;
        dY=(Yh-Yn.')/M; dW2=dY*A1.'; db2=sum(dY,2);
        dA1=W2.'*dY; dZ1=dA1.*(1-A1.^2); dW1=dZ1*Xt.'; db1=sum(dZ1,2);
        W1=W1-lr*dW1; b1=b1-lr*db1; W2=W2-lr*dW2; b2=b2-lr*db2;
    end
    model.type='ffnn'; model.W1=W1; model.b1=b1; model.W2=W2; model.b2=b2;
    model.xmu=xmu; model.xsig=xsig; model.ymu=ymu; model.ysig=ysig;
end

function F_cmd = nn_evaluate(model, pos_error, vel_error)
    if isempty(model), F_cmd=[0;0]; return; end
    ex=max(-0.4,min(0.4,pos_error(1))); ey=max(-0.4,min(0.4,pos_error(2)));
    vx=max(-0.3,min(0.3,vel_error(1))); vy=max(-0.3,min(0.3,vel_error(2)));
    Fx=nn_forward(model,[ex vx]); Fy=nn_forward(model,[ey vy]);
    F_cmd=[Fx;Fy];
    F_cmd(~isfinite(F_cmd))=0;
    F_cmd=max(-25,min(25,F_cmd));
end

function y = nn_forward(model, xrow)
    switch model.type
        case 'anfis'
            xrow(1)=max(model.in_range(1,1),min(model.in_range(1,2),xrow(1)));
            xrow(2)=max(model.in_range(2,1),min(model.in_range(2,2),xrow(2)));
            y=evalfis(model.fis,xrow);
            if ~isfinite(y), y=0; end
        case 'ffnn'
            xn=(xrow-model.xmu)./model.xsig;
            a1=tanh(model.W1*xn.'+model.b1); y=(model.W2*a1+model.b2)*model.ysig+model.ymu;
            if ~isfinite(y), y=0; end
        otherwise, y=0;
    end
end

function [X,Y] = build_nn_dataset(it2_fis, it2_scale, body_g, P, ref, vel)
    X=[]; Y=[];
    eg=linspace(-0.4,0.4,21); vg=linspace(-0.3,0.3,15);
    for b=1:numel(body_g)
        for ie=1:numel(eg)
            for iv=1:numel(vg)
                F=it2_evaluate(it2_fis,[eg(ie) 0],[vg(iv) 0],it2_scale(b));
                X=[X; eg(ie) vg(iv)]; Y=[Y; F(1)]; %#ok<AGROW>
            end
        end
    end
    conds={struct('mass_factor',1,'stiffness_factor',1), ...
           struct('mass_factor',1,'stiffness_factor',1), ...
           struct('mass_factor',1+P.mass_uncertainty,'stiffness_factor',1-P.stiffness_uncertainty)};
    inj=[false true false];
    for b=1:numel(body_g)
        Atmp.it2_fis=it2_fis; Atmp.it2_scale=it2_scale;
        for c=1:numel(conds)
            log=run_sim('IT2',b,Atmp,P,body_g(b),ref,vel,inj(c),conds{c});
            for k=1:5:P.N
                e_pos=log.error(k,:); ri=min(k,size(ref,1));
                e_vel=vel(ri,:)-log.vel_est(k,:);
                F=it2_evaluate(it2_fis,e_pos,e_vel,it2_scale(b));
                X=[X; e_pos(1) e_vel(1); e_pos(2) e_vel(2)]; Y=[Y; F(1); F(2)]; %#ok<AGROW>
            end
        end
    end
end

function F = it2_eval_wrap(asset,e,ev,s), F = it2_evaluate(asset,e,ev,s); end
function F = t3_eval_wrap(asset,e,ev,s),  F = t3_evaluate(asset,e,ev,s);  end

% ----- v3 robust output-scale tuner: 1-D grid scan + local refine -----
% returns the cost-minimising scale plus the coarse scan (for the figure).
% a grid search on a single scalar cannot oscillate/diverge the way the v2
% gradient descent did on the T3 surface, so IT2 and T3 both get their true
% optimum and the comparison is fair.
function [scale, scan] = tune_scale(evalfun, asset, ref_traj, P, g)
    lo=0.2; hi=3.0;
    coarse=lo:0.1:hi;
    cc=zeros(size(coarse));
    for ii=1:numel(coarse)
        cc(ii)=scale_cost(evalfun,asset,ref_traj,P,g,coarse(ii));
    end
    [~,bi]=min(cc); s_best=coarse(bi);
    % local refine around the coarse minimum
    fine=max(lo,s_best-0.15):0.02:min(hi,s_best+0.15);
    fc=zeros(size(fine));
    for ii=1:numel(fine)
        fc(ii)=scale_cost(evalfun,asset,ref_traj,P,g,fine(ii));
    end
    [~,fi]=min(fc); scale=fine(fi);
    scan.scales=coarse; scan.costs=cc;
    fprintf('.');
end

function cost = scale_cost(evalfun, asset, ref_traj, P, g, scale)
    Nr=size(ref_traj,1); pos=ref_traj(1,:); vel=[0 0];
    F_weight=[0,-P.mass_EE*g]; cost=0;
    L_actual=sqrt(sum((P.anchors-pos).^2,2)); alpha=min(P.dt*P.cable_compliance_bandwidth*2*pi,1.0);
    for k=1:Nr-1
        e_pos=ref_traj(k+1,:)-pos; e_vel=[0 0]-vel;
        F_ctrl=reshape(evalfun(asset,e_pos,e_vel,scale),1,[]);
        F_ext=F_weight+F_ctrl;
        [T_cmd,~]=tension_qp(pos,P.anchors,F_ext,P.T_min,P.T_max);
        [L_cmd,~,~]=ik_with_sag(pos,P.anchors,T_cmd,g,P.cable_mass_per_m);
        L_actual=L_actual+alpha*(L_cmd-L_actual); dL=L_cmd-L_actual;
        T_comp=max(P.T_min,min(P.T_max,T_cmd+P.cable_stiffness*dL-P.cable_damping*dL/P.dt));
        J=compute_jacobian(pos,P.anchors); F_net=(J'*T_comp)'+F_weight;
        acc=F_net/P.mass_EE; vel=vel+acc*P.dt; vel=max(-P.vel_max,min(P.vel_max,vel));
        pos=pos+vel*P.dt+0.5*acc*P.dt^2; pos=max(P.ws_min,min(P.ws_max,pos));
        cost=cost+sum(e_pos.^2);
    end
    if ~isfinite(cost), cost=1e9; end
end

% ----- generic simulation engine: dispatches to any of the four controllers -----
function log = run_sim(kind, b_idx, A, P, g, ref_traj, vel_traj, inject_dist, perturb)
    mass_actual=P.mass_EE*perturb.mass_factor;
    stiff_actual=P.cable_stiffness*perturb.stiffness_factor;
    F_weight=[0,-mass_actual*g];

    fo=[];
    if strcmpi(kind,'FOPID')
        sc=sqrt(g/9.81); pr=A.fopid_base; pr(1:3)=pr(1:3)*sc;
        fo=fopid_init(pr,P.dt,20.0,256);
    end

    pos=ref_traj(1,:); vel=[0 0];
    L_actual=sqrt(sum((P.anchors-pos).^2,2)); alpha=min(P.dt*P.cable_compliance_bandwidth*2*pi,1.0);
    N=P.N;
    log.t=P.t_vec; log.pos_ref=ref_traj;
    log.error=zeros(N,2); log.pos_est=zeros(N,2); log.vel_est=zeros(N,2);
    log.F_cmd=zeros(N,2);
    log.T_cables=zeros(N,P.num_cables); log.L_arc=zeros(N,P.num_cables);
    log.L_actual=zeros(N,P.num_cables); log.disturbance=zeros(N,2); log.stable=true;

    for k=1:N
        t=P.t_vec(k); ri=min(k,size(ref_traj,1));
        pos_des=ref_traj(ri,:); vel_des=vel_traj(ri,:);
        e_pos=pos_des-pos; e_vel=vel_des-vel;

        switch upper(kind)
            case 'IT2',   F_ctrl=it2_evaluate(A.it2_fis,e_pos,e_vel,A.it2_scale(b_idx));
            case 'T3',    F_ctrl=t3_evaluate(A.t3,e_pos,e_vel,A.t3_scale(b_idx));
            case 'FOPID', [F_ctrl,fo]=fopid_update(fo,pos_des,pos);
            case 'NN',    F_ctrl=nn_evaluate(A.nn_model,e_pos,e_vel);
            otherwise,    F_ctrl=[0;0];
        end
        F_ctrl=F_ctrl(:);
        F_ctrl(~isfinite(F_ctrl))=0;

        F_dist=[0;0];
        if inject_dist && t>=P.disturbance_time && t<=P.disturbance_time+P.disturbance_duration
            F_dist=[P.disturbance_magnitude;0];
        end
        A1=P.anchors(1,:); d1=norm(pos-A1);
        if d1>0.001, F_dist=F_dist+(P.tube_force_per_meter*d1*((A1-pos)/d1))'; end

        F_total_ext=F_weight'+F_dist; F_combined=F_total_ext(:)+F_ctrl;
        [T_cmd,~]=tension_qp(pos,P.anchors,F_combined,P.T_min,P.T_max);
        [L_cmd,~,~]=ik_with_sag(pos,P.anchors,T_cmd,g,P.cable_mass_per_m);
        L_actual=L_actual+alpha*(L_cmd-L_actual); dL=L_cmd-L_actual;
        T_comp=max(P.T_min,min(P.T_max,T_cmd+stiff_actual*dL-P.cable_damping*dL/P.dt));
        J=compute_jacobian(pos,P.anchors); F_net=(J'*T_comp)'+F_weight+F_dist';
        acc=F_net/mass_actual;
        vel=vel+acc*P.dt; vel=max(-1e4,min(1e4,vel));
        pos=pos+vel*P.dt+0.5*acc*P.dt^2; pos=max(P.ws_min,min(P.ws_max,pos));

        log.error(k,:)=e_pos; log.pos_est(k,:)=pos; log.vel_est(k,:)=vel;
        log.F_cmd(k,:)=F_ctrl'; log.T_cables(k,:)=T_cmd';
        log.L_arc(k,:)=L_cmd'; log.L_actual(k,:)=L_actual'; log.disturbance(k,:)=F_dist';

        if any(~isfinite(pos)) || any(~isfinite(vel))
            log.stable=false; log.error(k:end,:)=NaN; log.pos_est(k:end,:)=NaN; break;
        end
    end
    if max(sqrt(sum(log.error.^2,2))) > 2.0, log.stable=false; end
end

% ----- extended metrics -----
function m = compute_metrics_ext(log_step, log_dist, P)
    e_mag=sqrt(sum(log_step.error.^2,2))*1000;
    if all(~isfinite(e_mag))
        m=struct('rmse',NaN,'max_err',NaN,'sse',NaN,'overshoot',NaN,'settling',NaN, ...
                 'ISE',NaN,'IAE',NaN,'recovery',NaN,'tension_viol',NaN,'comp_lag',NaN, ...
                 'runtime',NaN,'stable',false); return;
    end
    e_mag(~isfinite(e_mag))=max(e_mag(isfinite(e_mag)));
    m.rmse=rms(e_mag); m.max_err=max(e_mag);
    m.sse=mean(e_mag(end-round(0.5/P.dt):end));
    ref_final=log_step.pos_ref(end,:); ref_dir=ref_final./norm(ref_final);
    err_dir=log_step.error*ref_dir';
    travel_mm=norm(log_step.pos_ref(end,:)-log_step.pos_ref(1,:))*1000;
    m.overshoot=max(0,-min(err_dir)*1000)/max(travel_mm,eps)*100;
    band=0.02*travel_mm; m.settling=P.t_total;
    for k=1:P.N, if all(e_mag(k:end)<band), m.settling=P.t_vec(k); break; end, end
    m.ISE=sum(e_mag.^2)*P.dt; m.IAE=sum(e_mag)*P.dt;
    ed=sqrt(sum(log_dist.error.^2,2))*1000; idx=find(P.t_vec>=P.disturbance_time,1);
    if ~isempty(idx) && idx<P.N
        ri=find(ed(idx:end)<5.0,1); if isempty(ri), m.recovery=NaN; else, m.recovery=ri*P.dt; end
    else, m.recovery=NaN;
    end
    m.tension_viol=sum(any(log_step.T_cables<P.T_min,2));
    dLa=log_step.L_arc-log_step.L_actual; m.comp_lag=mean(abs(dLa(:)))*1000;
    m.stable=log_step.stable && log_dist.stable;
end

% ----- springer figure helpers (copied from the lowgrav file) -----
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
