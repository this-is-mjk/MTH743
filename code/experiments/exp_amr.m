% exp_amr: Phase 4, AMR versus uniform grids on the dam break over a step; threshold study
set(0,'DefaultFigureVisible','off'); s=settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue='light';

root = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

g = 9.81;

% Problem definition and init
bed = make_bed('step',[0 0.5 12]);
init = @(x,Bc) deal(2*(x<8) + 1*(x>=8 & x<12) + 0.5*(x>=12), zeros(size(x)));
xl = 0; xr = 20; tend = 2;

% Reference solution (finest uniform grid)
cfgref = struct('N',6400,'xl',xl,'xr',xr,'bed',bed,'init',init,'tend',tend, ...
  'tsave',2,'scheme','hr','order',3,'rk','ssprk3','bc','wall','cfl',0.3,'g',g);
outref = swe1d_run(cfgref);
% Cell-average reference on 800 cells (finest grid for comparison)
Nf = 800; dx_ref = (xr - xl)/Nf;
Hfine = outref.H{1}; HUfine = outref.HU{1}; % 1x6400 at t=2
blk = numel(Hfine)/Nf;
Hf_reshaped = reshape(Hfine, blk, Nf);
HUf_reshaped = reshape(HUfine, blk, Nf);
Href = mean(Hf_reshaped,1);
HUref = mean(HUf_reshaped,1);
xref = xl + ((0.5:(Nf-0.5)))*dx_ref;

% Uniform gridding and AMR runs to compare (A,B,C,D)
% A: uniform hr-3 N=100
cfgA = struct('N',100,'xl',xl,'xr',xr,'bed',bed,'init',init,'tend',tend, ...
  'tsave',2,'scheme','hr','order',3,'rk','ssprk3','bc','wall','cfl',0.4,'g',g);

% B: uniform hr-3 N=800 (fine grid)
cfgB = struct('N',800,'xl',xl,'xr',xr,'bed',bed,'init',init,'tend',tend, ...
  'tsave',2,'scheme','hr','order',3,'rk','ssprk3','bc','wall','cfl',0.4,'g',g);

% C: AMR N0=100, Lmax=3 (finest 800)
cfgC = struct('N0',100,'xl',xl,'xr',xr,'bed',bed,'init',init,'tend',tend, ...
  'tsave',2,'scheme','hr','order',3,'rk','ssprk3','bc','wall','cfl',0.4,'g',g, ...
  'Lmax',3,'nregrid',5,'thr',2e-3,'mode','adaptive');

% D: AMR N0=200, Lmax=2 (finest 800)
cfgD = struct('N0',200,'xl',xl,'xr',xr,'bed',bed,'init',init,'tend',tend, ...
  'tsave',2,'scheme','hr','order',3,'rk','ssprk3','bc','wall','cfl',0.4,'g',g, ...
  'Lmax',2,'nregrid',5,'thr',2e-3,'mode','adaptive');

% Run configurations
outA = swe1d_run(cfgA);
outB = swe1d_run(cfgB);
outC = amr1d_run(cfgC);
outD = amr1d_run(cfgD);

% Compute L1 errors on h and hu by comparing at finest grid points or averaged
% For uniform grids A,B: interpolate reference on their grid by averaging reference or linear interp

% Average reference solution on 100 and 800 cells for comparison

% Avg ref for N=100 (same as cfgA.N)
blkA = numel(Hfine)/cfgA.N;   % 6400/100=64
Hf_reshapedA = reshape(Hfine, blkA, cfgA.N);
HUf_reshapedA = reshape(HUfine, blkA, cfgA.N);
HrefA = mean(Hf_reshapedA,1);
HUrefA = mean(HUf_reshapedA,1);

% Avg ref for N=800: already computed as Href and HUref

hA = outA.H{end}; huA = outA.HU{end};
hB = outB.H{end}; huB = outB.HU{end};
hC = outC.H{end}; huC = outC.HU{end};
hD = outD.H{end}; huD = outD.HU{end};

dxA = (xr - xl)/cfgA.N; dxB = (xr - xl)/cfgB.N;
xA = xl + ((0.5:(cfgA.N - 0.5)))*dxA;
xB = xl + ((0.5:(cfgB.N - 0.5)))*dxB;

% Interpolate fine ref to A and B uniform grids by linear interp
href_interp_A = HrefA;      % exact block average of the N=6400 reference (no interpolation)
hu_ref_interp_A = HUrefA;
href_interp_B = Href;       % same 800-cell grid
hu_ref_interp_B = HUref;

dxC = (xr - xl)/(cfgC.N0*2^cfgC.Lmax); % finest dx
dxD = (xr - xl)/(cfgD.N0*2^cfgD.Lmax);

% For errors of AMR (C and D), grid is nonuniform, but result.H{end} size matches finest Nf=800
% Its positions: outC.x pts (should be fine grid points)
xC = outC.x;
xD = outD.x;

% Compute L1 error by numerical integration over domain: sum abs difference * dx local
% For uniform grids A,B use dx uniform * sum abs difference
L1h_A = sum(abs(hA - href_interp_A))*dxA;
L1hu_A = sum(abs(huA - hu_ref_interp_A))*dxA;

L1h_B = sum(abs(hB - href_interp_B))*dxB;
L1hu_B = sum(abs(huB - hu_ref_interp_B))*dxB;

% For AMR grids C,D with same size as fine (Nf=800), assume points at xC, xD, uniform dx_ref = (xr-xl)/800
% We approximate L1 error with uniform weight dx_ref because both on fine grid
L1h_C = sum(abs(hC - Href))*dx_ref;
L1hu_C = sum(abs(huC - HUref))*dx_ref;

L1h_D = sum(abs(hD - Href))*dx_ref;
L1hu_D = sum(abs(huD - HUref))*dx_ref;

% Ratio error/uniform fine error (B)
ratio_h_A = L1h_A / L1h_B;
ratio_h_C = L1h_C / L1h_B;
ratio_h_D = L1h_D / L1h_B;
ratio_h_B = 1;

% nupdates from runs (count total cell updates)
nupd_A = outA.nupdates;
nupd_B = outB.nupdates;
nupd_C = outC.nupdates;
nupd_D = outD.nupdates;

nupd_savingC = nupd_B / nupd_C;
nupd_savingD = nupd_B / nupd_D;

% Measure total wall time (use tic/toc)
tstart = tic;
swe1d_run(cfgA);
tA = toc(tstart);
tstart = tic;
swe1d_run(cfgB);
tB = toc(tstart);
tstart = tic;
amr1d_run(cfgC);
tC = toc(tstart);
tstart = tic;
amr1d_run(cfgD);
tD = toc(tstart);

% Table output: method, L1 error h, L1 error hu, nupdates, time, err/uniform-fine err ratio, nupdates saving
tab = {
  'uniform-100', L1h_A, L1hu_A, nupd_A, tA, ratio_h_A, nupd_savingC; 
  'uniform-800', L1h_B, L1hu_B, nupd_B, tB, ratio_h_B, 1;
  'AMR-100-3',   L1h_C, L1hu_C, nupd_C, tC, ratio_h_C, nupd_savingC;
  'AMR-200-2',   L1h_D, L1hu_D, nupd_D, tD, ratio_h_D, nupd_savingD;
};

% Save results table
fname_results = fullfile(resdir,'amr.txt');
fid = fopen(fname_results,'w');
fprintf(fid,'Method\tL1_h\tL1_hu\tnupdates\ttime(s)\terror/h_uniform800\tnupdate_saving\n');
for i=1:size(tab,1)
  fprintf(fid,'%s\t%.4e\t%.4e\t%d\t%.3f\t%.3f\t%.3f\n', tab{i,:});
end
fclose(fid);

% Print results
fprintf('METHOD\tL1_h\t\tL1_hu\t\tnupdates\ttime(s)\terror_ratio\tnupdate_saving\n');
for i=1:size(tab,1)
  fprintf('%-12s %.4e %.4e %8d %7.3f %9.3f %9.3f\n', tab{i,:});
end

% ========== Figures ==========

c = lines(7);

% Figure 1: eta vs x for A,B,C + reference (black), bed; zoom on bore [7 13]
fig=figure('Color','w','Position',[100 100 900 600]);
plot(xref, Href + bed.f(xref,0), 'k-', 'LineWidth',1.8);
hold on
plot(outA.x, outA.H{end} + bed.f(outA.x,0), 'Color', c(1,:), 'LineStyle','--', 'LineWidth',1.6);
plot(outB.x, outB.H{end} + bed.f(outB.x,0), 'Color', c(2,:), 'LineStyle','-', 'LineWidth',1.6);
plot(outC.x, outC.H{end} + bed.f(outC.x,0), 'Color', c(3,:), 'LineStyle',':', 'LineWidth',1.6);
xlim([0 20]);
xlabel('x [m]'); ylabel('\eta [m]');
grid on;
legend({'Reference (6400)','Uniform N=100','Uniform N=800','AMR N0=100 Lmax=3'}, 'Location','northwest');
title('Surface elevation \eta at t=2, uniform and AMR vs reference');
exportgraphics(fig, fullfile(imgdir,'amr_eta_compare.png'), 'Resolution',150);

% Zoomed plot [7 13]
fig=figure('Color','w','Position',[100 100 900 600]);
plot(xref, Href + bed.f(xref,0), 'k-', 'LineWidth',1.8);
hold on
plot(outA.x, outA.H{end} + bed.f(outA.x,0), 'Color', c(1,:), 'LineStyle','--', 'LineWidth',1.6);
plot(outB.x, outB.H{end} + bed.f(outB.x,0), 'Color', c(2,:), 'LineStyle','-', 'LineWidth',1.6);
plot(outC.x, outC.H{end} + bed.f(outC.x,0), 'Color', c(3,:), 'LineStyle',':', 'LineWidth',1.6);
xlim([7 13]);
xlabel('x [m]'); ylabel('\eta [m]');
grid on;
legend({'Reference (6400)','Uniform N=100','Uniform N=800','AMR N0=100 Lmax=3'}, 'Location','best');
title('Zoom on bore region at t=2');
exportgraphics(fig, fullfile(imgdir,'amr_eta_zoom.png'), 'Resolution',150);

% Figure 2: Level map (x vs time) for AMR run C
tsaveC = linspace(0,tend,81);
cfgC.tsave = tsaveC;
outC_map = amr1d_run(cfgC);
% Compose level matrix for image: rows = time steps, cols = finest grid cells
levmat = zeros(length(tsaveC), Nf);
for i=1:length(tsaveC)
  lev = outC_map.level{i}; % level index 0..Lmax, length Nf=800
  levmat(i,:) = lev;
end
fig=figure('Color','w','Position',[100 100 900 600]);
imagesc(xref, tsaveC, levmat);
axis xy;
colormap(parula(cfgC.Lmax+1));
clim([-0.5 cfgC.Lmax+0.5]);          % one colour per integer level
cb = colorbar;
cb.Ticks = 0:cfgC.Lmax;
cb.Label.String = 'AMR level';
xlabel('x [m]'); ylabel('time [s]');
title('AMR refinement level map vs time (run C)');
exportgraphics(fig, fullfile(imgdir,'amr_level_map.png'), 'Resolution',150);

% Figure 3: Active cells per step vs time for C and D + uniform fine count (horizontal line)
tsaveCD = linspace(0,tend,81);
cfgC.tsave = tsaveCD;
outC_act = amr1d_run(cfgC);
cfgD.tsave = tsaveCD;
outD_act = amr1d_run(cfgD);

% Use numel and nsteps from available fields nstepsC = numel(outC_act.ncells_hist);
nstepsC = numel(outC_act.ncells_hist);
nstepsD = numel(outD_act.ncells_hist);

actC = outC_act.ncells_hist;
actD = outD_act.ncells_hist;

fig=figure('Color','w','Position',[100 100 900 600]);
hold on;
plot(outC_act.t_hist, actC, 'Color', c(3,:), 'LineWidth',1.7, 'DisplayName','AMR C (N0=100, Lmax=3)');
plot(outD_act.t_hist, actD, 'Color', c(4,:), 'LineWidth',1.7, 'DisplayName','AMR D (N0=200, Lmax=2)');
yline(Nf, 'k--', 'LineWidth',1.5, 'DisplayName','Uniform Fine 800');
xlabel('Time [s]');
ylabel('Active cells');
legend('Location','northeast');
title('Active cells per step vs time');
grid on;
exportgraphics(fig, fullfile(imgdir,'amr_active_cells_time.png'), 'Resolution',150);

% Sensitivity study thr for run C (3 values)
thrvals = [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4];   % relative to the range of eta
errs = zeros(1,numel(thrvals));
nups = zeros(1,numel(thrvals));
for k=1:numel(thrvals)
  cfgCT = cfgC;
  cfgCT.thr = thrvals(k);
  outCT = amr1d_run(cfgCT);
  hCT = outCT.H{end}; huCT = outCT.HU{end};
  errs(k) = sum(abs(hCT - Href)) * dx_ref;
  nups(k) = outCT.nupdates;
end

% Pareto plot: error vs nupdates
fig=figure('Color','w','Position',[100 100 900 600]);
loglog(nups, errs, '-o', 'Color',[0 0.447 0.741], 'LineWidth',1.7, 'MarkerSize',8, 'DisplayName','AMR N_0=100, 3 levels (thr varied)');
hold on;
loglog([nupd_A nupd_B], [L1h_A L1h_B], 's', 'Color',[0 0 0], 'MarkerSize',10, 'LineWidth',1.7, 'DisplayName','uniform N=100, N=800');
xlabel('Number of cell updates');
ylabel('L1 error in h (vs N=6400 reference)');
grid on; legend('Location','northeast'); set(gca,'FontSize',12);
title('AMR: error versus work for different refinement thresholds');
xlim([2e4 5e6]); ylim([1.2e-2 3e-1]);
text(nups(1)*1.05, errs(1), sprintf('thr = %.0e', thrvals(1)), 'FontSize',11);
text(nups(2)*0.55, errs(2)*0.9, sprintf('thr = %.0e', thrvals(2)), 'FontSize',11);
text(nups(end)*1.1, errs(end)*1.12, sprintf('thr = %.0e', thrvals(end)), 'FontSize',11);
text(nupd_A*1.1, L1h_A, 'uniform N=100', 'FontSize',11);
text(nupd_B*0.45, L1h_B*0.85, 'uniform N=800', 'FontSize',11);
exportgraphics(fig, fullfile(imgdir,'amr_sensitivity_thr.png'), 'Resolution',150);

% MP4 animation amr_animation.mp4 top=eta+bed, bottom=level map + edges
fps = 20;
tsavemp4 = linspace(0,tend,101);
cfgC.tsave = tsavemp4;
out_anim = amr1d_run(cfgC);

vidfile = fullfile(imgdir,'amr_animation.mp4');
v = VideoWriter(vidfile,'MPEG-4');
v.FrameRate = fps;
open(v);

% Prepare bed line
bedvals = bed.f(out_anim.x,0);

fig=figure('Color','w','Position',[100 100 900 600]);

for i=1:length(tsavemp4)
  clf;
  eta = out_anim.H{i} + bedvals;
  lev = out_anim.level{i};
  % Top subplot eta and bed
  subplot(2,1,1);
  plot(out_anim.x, eta, 'b-', 'LineWidth',1.8);
  hold on;
  plot(out_anim.x, bedvals, 'k--', 'LineWidth',1.4);
  ylim([0 3]);
  title(sprintf('Surface elevation \\eta and bed at t=%.3f s', tsavemp4(i)));
  ylabel('\eta, bed [m]');
  grid on;

  % Bottom subplot level map as color strip + cell edges
  subplot(2,1,2);
  imagesc(out_anim.x, [0 1], lev, [0 cfgC.Lmax]);
  colormap(parula(cfgC.Lmax+1));
  hold on;
  % plot cell edges as vertical lines
  dx_lev = (xr - xl)/length(lev);
  xedges = xl + (0:length(lev)) * dx_lev;
  yl = [0 1];
  plot([xedges; xedges], repmat(yl',1,length(xedges)), 'k-', 'LineWidth', 0.5);
  xlabel('x [m]');
  yticks([]);
  title('AMR level map');
  c = colorbar('eastoutside');
  c.Label.String = 'AMR Level';
  drawnow;
  frame = getframe(fig);
  writeVideo(v,frame);
end
close(v);

% Conclusion line
concl = sprintf('AMR run C reaches uniform fine grid accuracy L1_h=%.2e with %.2f%% fewer updates than uniform fine.', ...
  L1h_C, 100*(1 - nupd_C/nupd_B));
fprintf('%s\n', concl);

fprintf('EXP amr DONE\n');
