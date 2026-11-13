function exp_rk
% RUN stability, accuracy, efficiency study for SWE 1D time integrators
% with cosbump bed, wall bc, hr scheme order 3, N=100, tend=0.5
% Constraint: <4 min runtime, no toolboxes, vectorized

set(0,'DefaultFigureVisible','off');
s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue='light';

root = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

% Config base problem parameters
N = 100; xl=0; xr=1;
bed = make_bed('cosbump',[0.5 0.5 0.1]);
bc = 'wall'; scheme = 'hr'; order = 3;
tend = 0.5; % short time, enough for instabilities
init = @(x,Bc) deal(max(0, 1 + 0.2*(x>0.1 & x<0.2) - Bc), zeros(size(x)));

% Integrator list and info: stages, order, SSP, storage regs
rk_list = {'euler','ssprk2','ssprk3','rk4'};
stages = [1 2 3 4];
formal_order = [1 2 3 4];
SSP = [true true true false];
storage = [1 2 2 3];

% Reference solution: ssprk3 cfl=0.05 on same N, store h,u max velocity
cfg_ref = struct();
cfg_ref.N = N; cfg_ref.xl = xl; cfg_ref.xr = xr;
cfg_ref.bed = bed;
cfg_ref.bc = bc;
cfg_ref.scheme = scheme;
cfg_ref.order = order;
cfg_ref.rk = 'ssprk3';
cfg_ref.cfl = 0.05;
cfg_ref.tend = tend;
cfg_ref.init = init;
cfg_ref.tsave = tend;
cfg_ref.monitor = true;
% Run reference
fprintf('Running reference solution (ssprk3, CFL=0.05)... ');
out_ref = swe1d_run(cfg_ref);
fprintf('done.\n');
% Check max|u| at tend > 0.1 to confirm pulse moves
u_ref = out_ref.HU{end} ./ out_ref.H{end};
u_ref(isnan(u_ref) | isinf(u_ref)) = 0; % guard unlikely artifacts
max_vel_ref = max(abs(u_ref));
fprintf('Reference max |u| at t=%.2f = %.5f\n', tend, max_vel_ref);
if max_vel_ref<=0.1
    warning('Reference max velocity <=0.1 might indicate no pulse movement');
end

% Precompute dx for error calc
dx = (xr - xl)/N;

% --- 1. Stability: find max stable CFL in [0.3,4] (start at 0.3) ---
% Method stable if swe1d_run finishes with finite H, max(H+B)<1.5, min(H)>0
cfl_start = 0.3;
cfl_maxscan = 4;
cfl_factor = 1.15;
bisection_steps = 6;
cfl_low_bound = 0.02;

stableCFL = nan(size(rk_list));
fprintf('--- Stability scan ---\n');
for i=1:numel(rk_list)
    rk = rk_list{i};
    % Scan CFL upwards from cfl_start
    cfl_try = cfl_start;
    % Keep track of last stable and first unstable CFL for bisection
    cfl_stable = NaN;
    cfl_unstable = NaN;
    stablefound = false;
    % Special: test if unstable even at cfl=0.02
    unstableAtLow = ~run_stability_test(N, xl, xr, bed, bc, scheme, order, rk, 0.02, tend, init);
    if unstableAtLow
        fprintf('Integrator %s unstable for all CFL\n', rk);
        stableCFL(i) = NaN; % Mark no stable CFL
        continue;
    end
    % Scan increasing CFL
    while cfl_try <= cfl_maxscan
        ok = run_stability_test(N, xl, xr, bed, bc, scheme, order, rk, cfl_try, tend, init);
        if ok
            cfl_stable = cfl_try;
            cfl_try = cfl_try * cfl_factor;
        else
            cfl_unstable = cfl_try;
            stablefound = true;
            break;
        end
    end
    if isnan(cfl_unstable)
        % Up to max scan all stable
        fprintf('Integrator %s limit above scan range\n', rk);
        stableCFL(i) = Inf;
    else
        % Bisection refinement between cfl_stable and cfl_unstable
        a = cfl_stable;
        b = cfl_unstable;
        for bs=1:bisection_steps
            m = (a+b)/2;
            okm = run_stability_test(N, xl, xr, bed, bc, scheme, order, rk, m, tend, init);
            if okm
                a = m;
            else
                b = m;
            end
        end
        stableCFL(i) = a;
        fprintf('Integrator %s max stable CFL = %.6f\n', rk, a);
    end
end

% --- 2. Accuracy: L1 error in h vs ref at cfl=0.4 and 0.8*max_CFL ---
% Only at stable CFLs
fprintf('--- Accuracy tests ---\n');
cfls_acc = nan(2,numel(rk_list)); % rows: normal (0.4) and scaled (0.8 max)
errors = nan(2,numel(rk_list));
for i=1:numel(rk_list)
    rk = rk_list{i};
    maxCFL = stableCFL(i);
    if isnan(maxCFL) || maxCFL<0.4
        % Too unstable, skip
        fprintf('Skipping accuracy for %s (no stable CFL >=0.4)\n', rk);
        continue;
    end
    % 1) test cfl=0.4
    c1 = 0.4;
    cfls_acc(1,i) = c1;
    err1 = run_error(N, xl, xr, bed, bc, scheme, order, rk, c1, tend, init, out_ref);
    errors(1,i) = err1;
    % 2) test cfl=0.8*max stable CFL
    c2 = 0.8*maxCFL;
    cfls_acc(2,i) = c2;
    err2 = run_error(N, xl, xr, bed, bc, scheme, order, rk, c2, tend, init, out_ref);
    errors(2,i) = err2;
    fprintf('Integrator %s error at CFL=%.3f: %.5e, at CFL=%.3f: %.5e\n', rk, c1, err1, c2, err2);
end

% --- 3. Efficiency: error vs nupdates (cost) curve ---
% CFL values: 0.3:0.15: max stable (max 6 values per method)
fprintf('--- Efficiency runs ---\n');
effData = cell(numel(rk_list),2); % {errors, nupdates} per method
for i=1:numel(rk_list)
    rk = rk_list{i};
    maxCFL = stableCFL(i);
    if isnan(maxCFL) || maxCFL < 0.3
        fprintf('Skipping efficiency for %s (no stable CFL >=0.3)\n', rk);
        effData{i,1} = [];
        effData{i,2} = [];
        continue;
    end
    % CFL values from 0.3 by 0.15 up to maxCFL, max 6 points
    CFLs = 0.3 : 0.15 : maxCFL;
    if numel(CFLs) > 6
        CFLs = CFLs(1:6);
    end
    errs = nan(size(CFLs));
    nups = nan(size(CFLs));
    for j=1:numel(CFLs)
        cflj = CFLs(j);
        [errj,nupj] = run_error_withcost(N, xl, xr, bed, bc, scheme, order, rk, cflj, tend, init, out_ref);
        errs(j) = errj;
        nups(j) = nupj;
        fprintf('Eff %s CFL=%.3f err=%.3e nupdates=%d\n', rk, cflj, errj, nupj);
    end
    effData{i,1} = errs;
    effData{i,2} = nups;
end

% --- 4. Create table with data ---
% Columns: integrator, stages, order, SSP, storage, max stable CFL (or Inf), err@0.4, err@0.8*max
fname_txt = fullfile(resdir,'exp_rk.txt');
fid = fopen(fname_txt,'w');
fprintf(fid, 'Integrator     Stages   Order   SSP   Storage  MaxCFL         Error@0.4    Error@0.8*max\n');
fprintf(fid, '---------------------------------------------------------------------------------------\n');
for i=1:numel(rk_list)
    rk = rk_list{i};
    sspstr = 'yes';
    if ~SSP(i), sspstr='no'; end
    maxCFLval = stableCFL(i);
    if isinf(maxCFLval), maxCFLstr = 'above scan'; else maxCFLstr=sprintf('%.6f',maxCFLval); end
    if isnan(errors(1,i)), e1s='-'; else e1s=sprintf('%.3e',errors(1,i)); end
    if isnan(errors(2,i)), e2s='-'; else e2s=sprintf('%.3e',errors(2,i)); end
    fprintf(fid, '%-13s %6d %7d %5s %8d %13s %14s %15s\n', ...
        rk, stages(i), formal_order(i), sspstr, storage(i), maxCFLstr, e1s, e2s);
end
fclose(fid);

% Also print to command window
type(fname_txt);

% --- 5. Plot efficiency: error vs nupdates loglog ---
colors = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556];   % no yellow (unreadable on white)
mk = {'o','s','^','d'};
fig = figure('Color','w','Position',[100 100 900 600]);
hold on; grid on; box on;
for i=1:numel(rk_list)
    e = effData{i,1};
    nup = effData{i,2};
    if isempty(e), continue; end
    plot(nup,e,mk{i},'MarkerSize',8,'LineWidth',1.8,'Color',colors(i,:),'MarkerFaceColor',colors(i,:));
end
xlabel('nupdates (cost)');
ylabel('L1 error in h');
set(gca,'XScale','log','YScale','log','FontSize',13);
legend(rk_list,'Location','northeast','Interpreter','none');
title('Error versus cost of the time integrators (hr order 3, N=100): lower left is better');
exportgraphics(fig, fullfile(imgdir,'rk_efficiency.png'), 'Resolution',150);

% --- Conclusion ---
% Least error per cost from final point (min error/nupdates)
min_ratio = Inf;
best_idx = 0;
for i=1:numel(rk_list)
    e = effData{i,1};
    nup = effData{i,2};
    if isempty(e), continue; end
    ratio = e(end)/nup(end);
    if ratio < min_ratio
        min_ratio = ratio;
        best_idx = i;
    end
end
if best_idx>0
    concl = sprintf('Integrator %s gives the least error per cost among tested methods.', rk_list{best_idx});
else
    concl = 'No integrator produced stable accurate solutions.';
end
fprintf('\n%s\n', concl);

fprintf('EXP rk DONE\n');

end

% === Helpers ===

function ok = run_stability_test(N, xl, xr, bed, bc, scheme, order, rk, cfl, tend, init)
% Runs problem with given rk and cfl, returns true if stable per criteria
cfg = struct();
cfg.N = N; cfg.xl = xl; cfg.xr = xr;
cfg.bed = bed; cfg.bc = bc; cfg.scheme = scheme; cfg.order = order;
cfg.rk = rk; cfg.cfl = cfl; cfg.tend = tend; cfg.init = init;
cfg.tsave = tend;
try
    out = swe1d_run(cfg);
catch
    ok = false; return;
end
H = out.H{end};
B = out.B;
ok = all(isfinite(H)) && max(H+B) < 1.5 && min(H) > 0;
end

function err = run_error(N, xl, xr, bed, bc, scheme, order, rk, cfl, tend, init, out_ref)
% Run and compute L1 error vs out_ref at tend
cfg = struct();
cfg.N = N; cfg.xl = xl; cfg.xr = xr;
cfg.bed = bed; cfg.bc = bc; cfg.scheme = scheme; cfg.order = order;
cfg.rk = rk; cfg.cfl = cfl; cfg.tend = tend; cfg.init = init;
cfg.tsave = tend;
out = swe1d_run(cfg);
hnum = out.H{end};
href = out_ref.H{end};
dx = (xr - xl)/N;
err = sum(abs(hnum - href))*dx;
end

function [err,nupdates] = run_error_withcost(N, xl, xr, bed, bc, scheme, order, rk, cfl, tend, init, out_ref)
% Run, return L1 error and nupdates count
cfg = struct();
cfg.N = N; cfg.xl = xl; cfg.xr = xr;
cfg.bed = bed; cfg.bc = bc; cfg.scheme = scheme; cfg.order = order;
cfg.rk = rk; cfg.cfl = cfl; cfg.tend = tend; cfg.init = init;
cfg.tsave = tend;
out = swe1d_run(cfg);
hnum = out.H{end};
href = out_ref.H{end};
dx = (xr - xl)/N;
err = sum(abs(hnum - href))*dx;
nupdates = out.nupdates;
end
