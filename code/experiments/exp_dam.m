% exp_dam: dam break over a step, hr-1 and hr-3 against a fine reference; animation
set(0,'DefaultFigureVisible','off'); s=settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue='light';

root = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures'); 
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

%% Problem setup
g = 9.81;
xl=0; xr=20;

% Bed: step with B0=0, B1=0.5 at xs=12
bed = make_bed('step',[0 0.5 12]);

% Initial water surface elevation eta0 stratified:
% eta = 2 for x<8 (B=0), 
% eta = 1 for 8<=x<12 (B=0), 
% eta = 1 for x>=12 (B=0.5) so h=eta-B=0.5 for x>=12 submerged by 0.5

init_fun = @(x,Bc) dam_init(x,Bc);

% Times to save for Figure C
tsave_fine = linspace(0,2,101); % for animation
tsave_coarse = [0.5 1 1.5 2];   % for x-t plot and output

% Common domain and times
tend = 2;

% Methods configurations: name, N, order, rk, cfl
configs = {
  struct('name','hr-1','N',400,'order',1,'rk','euler','cfl',0.5);
  struct('name','hr-3-400','N',400,'order',3,'rk','ssprk3','cfl',0.4);
  struct('name','hr-3-100','N',100,'order',3,'rk','ssprk3','cfl',0.4);
  struct('name','ref','N',6400,'order',3,'rk','ssprk3','cfl',0.3);
};

% Storage results
results = struct();

%% Run REFERENCE FIRST
for i=1:length(configs)
  cfg = configs{i};
  if strcmp(cfg.name,'ref')
    cfg.xl = xl; cfg.xr = xr; cfg.tend = tend;
    cfg.tsave = tsave_coarse;
    cfg.bed = bed;
    cfg.init = init_fun;
    cfg.scheme = 'hr';
    cfg.bc = 'wall';
    cfg.g = g;
    cfg.monitor = true;
    cfg.cfl = cfg.cfl;
    cfg.order = cfg.order;
    cfg.rk = cfg.rk;
    cfg.N = cfg.N;

    fprintf('Running %s N=%d order=%d rk=%s cfl=%.2g\n', cfg.name, cfg.N, cfg.order, cfg.rk, cfg.cfl);
    tic;
    out = swe1d_run(cfg);
    walltime = toc;

    % For reference run, zero errors
    err_h = 0;
    err_hu = 0;

    % Save results
    results.ref = struct('cfg',cfg,'out',out);

    results.ref.err_h = err_h;
    results.ref.err_hu = err_hu;
    results.ref.walltime = walltime;
    results.ref.nupdates = out.nupdates;

    % ref errors per cost NaN
    results.ref.err_per_cost_h = NaN;
    results.ref.err_per_cost_hu = NaN;

    % Compute initial mass from initial condition
    N = cfg.N;
    dx = (xr - xl)/N;
    [h0,u0] = init_fun(out.x, out.B);
    mass_init = sum(h0)*dx;
    results.ref.mass_init = mass_init;

    % Compute final mass from last snapshot H{end}
    mass_final = sum(out.H{end})*dx;
    results.ref.mass_final = mass_final;

    % Relative drift
    rel_drift = abs(mass_final - mass_init) / mass_init;
    results.ref.rel_drift = rel_drift;

    break;
  end
end

%% Then run others
for i=1:length(configs)
  cfg = configs{i};
  if strcmp(cfg.name, 'ref')
    continue
  end

  cfg.xl = xl; cfg.xr = xr; cfg.tend = tend;
  cfg.tsave = tsave_coarse;
  cfg.bed = bed;
  cfg.init = init_fun;
  cfg.scheme = 'hr';
  cfg.bc = 'wall';
  cfg.g = g;
  cfg.monitor = true;
  cfg.cfl = cfg.cfl;
  cfg.order = cfg.order;
  cfg.rk = cfg.rk;
  cfg.N = cfg.N;

  fprintf('Running %s N=%d order=%d rk=%s cfl=%.2g\n', cfg.name, cfg.N, cfg.order, cfg.rk, cfg.cfl);
  tic;
  out = swe1d_run(cfg);
  walltime = toc;

  % Compute L1 errors in h and hu w.r.t. ref (restricted to coarse grids)
  % Restrict ref solution by averaging
  refdata = results.ref;
  ratio = refdata.cfg.N / cfg.N;
  dx_ref = (xr - xl) / refdata.cfg.N;
  dx_coarse = (xr - xl) / cfg.N;
  % Average ref to coarse grid for error
  H_ref = refdata.out.H{end};
  HU_ref = refdata.out.HU{end};
  H_ref_coarse = zeros(1,cfg.N);
  HU_ref_coarse = zeros(1,cfg.N);
  for k=1:cfg.N
    ind_ref = (1:ratio) + (k-1)*ratio;
    H_ref_coarse(k) = mean(H_ref(ind_ref));
    HU_ref_coarse(k) = mean(HU_ref(ind_ref));
  end
  dx_c = dx_coarse;
  err_h = sum(abs(out.H{end} - H_ref_coarse))*dx_c;
  err_hu = sum(abs(out.HU{end} - HU_ref_coarse))*dx_c;

  % Save results
  fldname = strrep(cfg.name,'-','_');
  results.(fldname).out = out;
  results.(fldname).cfg = cfg;
  results.(fldname).err_h = err_h;
  results.(fldname).err_hu = err_hu;
  results.(fldname).walltime = walltime;
  results.(fldname).nupdates = out.nupdates;

  % Error per cost
  if err_h>0
    results.(fldname).err_per_cost_h = err_h / walltime;
  else
    results.(fldname).err_per_cost_h = NaN; % just in case
  end
  if err_hu>0
    results.(fldname).err_per_cost_hu = err_hu / walltime;
  else
    results.(fldname).err_per_cost_hu = NaN;
  end

  % Compute initial mass from initial condition
  N = cfg.N;
  dx = (xr - xl) / N;
  [h0,u0] = init_fun(out.x, out.B);
  mass_init = sum(h0) * dx;
  results.(fldname).mass_init = mass_init;

  % Compute final mass from last snapshot H{end}
  mass_final = sum(out.H{end}) * dx;
  results.(fldname).mass_final = mass_final;

  % Relative drift
  rel_drift = abs(mass_final - mass_init) / mass_init;
  results.(fldname).rel_drift = rel_drift;
end

%% Figure A: eta=h+B at t=2 with zoom inset x in [9,15]
figA=figure('Color','w','Position',[100 100 900 600]); hold on;
c=lines(7);

% Plot bed filled grey
xA = results.hr_3_400.out.x;
B = results.hr_3_400.out.B;
fill([xA fliplr(xA)],[B zeros(size(B))],'k','FaceAlpha',0.15,'EdgeColor','none','DisplayName','bed');

% Plot eta=h+B
plot(results.hr_1.out.x, results.hr_1.out.H{end}+B, 'Color', c(1,:), 'LineWidth', 1.8, 'LineStyle','-','DisplayName','hr-1 N=400')
plot(results.hr_3_400.out.x, results.hr_3_400.out.H{end}+B, 'Color', c(2,:), 'LineWidth', 1.8, 'LineStyle', '--','DisplayName','hr-3 N=400')
plot(results.hr_3_100.out.x, results.hr_3_100.out.H{end}+bed.f(results.hr_3_100.out.x,0), 'Color', c(3,:), 'LineWidth', 1.8, 'LineStyle', ':','DisplayName','hr-3 N=100')
plot(results.ref.out.x, results.ref.out.H{end}+results.ref.out.B, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '-', 'DisplayName','Reference N=6400')

xlabel('x [m]'), ylabel('\eta = h + B [m]')
title('Figure A: Water surface (\eta) at t=2; bed filled grey')
legend('Location','NorthWest')
xlim([0 20])
ylim([0 2.1])
grid on

% Inset zoom x in [9,15]
axes('Position',[0.16 0.16 0.34 0.30]); hold on     % empty lower-left corner of the main plot
fill([xA fliplr(xA)],[B zeros(size(B))],'k','FaceAlpha',0.15,'EdgeColor','none');
plot(results.hr_1.out.x, results.hr_1.out.H{end}+B, 'Color', c(1,:), 'LineWidth', 1.5, 'LineStyle','-')
plot(results.hr_3_400.out.x, results.hr_3_400.out.H{end}+B, 'Color', c(2,:), 'LineWidth', 1.5, 'LineStyle','--')
plot(results.hr_3_100.out.x, results.hr_3_100.out.H{end}+bed.f(results.hr_3_100.out.x,0), 'Color', c(3,:), 'LineWidth', 1.5, 'LineStyle', ':')
plot(results.ref.out.x, results.ref.out.H{end}+results.ref.out.B, 'Color', [0 0 0], 'LineWidth', 1.5, 'LineStyle','-')
xlim([9 15])
ylim([0 1.75])      % the surface is at 1.43-1.57 here: the old limit 1.4 cut every curve
grid on
box on

exportgraphics(figA, fullfile(imgdir,'dam_eta_t2.png'),'Resolution',150);

%% Figure B: u(x) at t=2 for all methods
figB=figure('Color','w','Position',[100 100 900 600]); hold on;
plot_u = @(out,B) safe_u(out.H{end}, out.HU{end});
plot(results.hr_1.out.x, plot_u(results.hr_1.out,results.hr_1.out.B), 'Color', c(1,:), 'LineWidth', 1.8, 'LineStyle','-','DisplayName','hr-1 N=400')
plot(results.hr_3_400.out.x, plot_u(results.hr_3_400.out,results.hr_3_400.out.B), 'Color', c(2,:), 'LineWidth', 1.8, 'LineStyle','--','DisplayName','hr-3 N=400')
plot(results.hr_3_100.out.x, plot_u(results.hr_3_100.out,results.hr_3_100.out.B), 'Color', c(3,:), 'LineWidth', 1.8, 'LineStyle',':','DisplayName','hr-3 N=100')
plot(results.ref.out.x, plot_u(results.ref.out,results.ref.out.B), 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle','-','DisplayName','Reference N=6400')

xlabel('x [m]'), ylabel('u [m/s]')
title('Figure B: Velocity u at t=2 for all methods')
legend('Location','NorthEast')
xlim([0 20])
grid on

exportgraphics(figB, fullfile(imgdir,'dam_u_t2.png'),'Resolution',150);

%% Figure C: x-t style stacked eta snapshots for hr-3 N=400 from 4 saved times
figC=figure('Color','w','Position',[100 100 900 600]);
hold on;
offset = 0;
for k=1:length(tsave_coarse)
  xC = results.hr_3_400.out.x;
  B = results.hr_3_400.out.B;
  eta = results.hr_3_400.out.H{k} + B;
  plot(xC, eta + offset, 'LineWidth', 1.8, 'DisplayName', sprintf('t=%.1f', tsave_coarse(k)));
  offset = offset + 1.2; % stack plots offset
end
xlabel('x [m]')
ylabel('Stacked \eta(t) snapshots [m]')
title('Figure C: Stacked water surface \eta(x) at saved times for hr-3 N=400')
legend('Location','NorthWest')
xlim([0 20])
grid on

exportgraphics(figC, fullfile(imgdir,'dam_eta_xt_stacked.png'),'Resolution',150);

%% Table: method, N, L1 err h, L1 err hu, nupdates, wall time, error per cost (h only)
tabfile = fullfile(resdir,'dam.txt');
fid = fopen(tabfile,'w');
fprintf(fid, 'Method\tN\tL1 err h\tL1 err hu\tnupdates\twall time [s]\terr per cost h\n');
fprintf('Method    N    L1 err h    L1 err hu    nupdates    wall time [s]    err per cost h\n');
fields_order = {'hr_1','hr_3_400','hr_3_100','ref'};
for k=1:length(fields_order)
  f = fields_order{k};
  r = results.(f);
  fprintf(fid, '%s\t%d\t%.4e\t%.4e\t%d\t%.3f\t%.4e\n', strrep(f,'_','-'), r.cfg.N, r.err_h, r.err_hu, r.nupdates, r.walltime, r.err_per_cost_h);
  fprintf('%-9s %4d %10.4e %10.4e %10d %12.3f %12.4e\n', strrep(f,'_','-'), r.cfg.N, r.err_h, r.err_hu, r.nupdates, r.walltime, r.err_per_cost_h);
end
fclose(fid);

%% Check conservation: relative total mass drift for all runs
% Initial mass sum(h0*dx)
% Adjusted calculation using initial condition function to compute mass_init
fprintf('Mass conservation relative drift:\n');
massfile = fullfile(resdir,'dam_mass.txt');
fidm = fopen(massfile,'w');
fprintf(fidm,'Method\tN\tMass0\tMass(tend)\tRelative drift\n');
for k=1:length(fields_order)
  f = fields_order{k};
  r = results.(f);

  % mass_init = mass0([r.H0; r.HU0], r.cfg); % Removed old call with non-existent fields
  mass_init = r.mass_init;
  mass_final = r.mass_final;
  rel_drift = r.rel_drift;

  fprintf(fidm, '%s\t%d\t%.6f\t%.6f\t%.3e\n', strrep(f,'_','-'), r.cfg.N, mass_init, mass_final, rel_drift);
  fprintf('%-9s %4d mass0=%.6f massf=%.6f rel drift=%.3e\n', strrep(f,'_','-'), r.cfg.N, mass_init, mass_final, rel_drift);
end
fclose(fidm);

%% MP4 Animation for hr-3 N=400, tsave=linspace(0,2,101), 20 fps
% We redo this run with finer outputs just for animation
anim_cfg = results.hr_3_400.cfg;
anim_cfg.tsave = tsave_fine;
anim_out = swe1d_run(anim_cfg);
fig_anim = figure('Color','w','Position',[100 100 900 700]);
nframes = length(tsave_fine);
filename = fullfile(root,'outputs','figures','dam_animation.mp4');
v = VideoWriter(filename,'MPEG-4');
v.FrameRate=20;
open(v);

for k=1:nframes
  clf;
  % Top subplot: eta and bed
  subplot(2,1,1);
  hold on;
  B = anim_out.B;
  x = anim_out.x;
  eta = anim_out.H{k}+B;
  plot(x,B,'Color',[0.5 0.5 0.5],'LineWidth',2,'DisplayName','Bed B')
  plot(x,eta,'b-','LineWidth',2,'DisplayName','Water surface \eta')
  xlabel('x [m]')
  ylabel('\eta, B [m]')
  title(sprintf('Water surface and bed at t=%.2f s', tsave_fine(k)))
  xlim([0 20])
  ylim([0 2.5])
  legend('Location','NorthWest')
  grid on
  
  % Bottom subplot: velocity u
  subplot(2,1,2);
  h = anim_out.H{k};
  hu = anim_out.HU{k};
  u = zeros(size(h));
  wet = h > 1e-10;
  u(wet) = hu(wet)./h(wet);
  plot(x,u,'r-','LineWidth',2)
  xlabel('x [m]')
  ylabel('u [m/s]')
  title('Velocity u')
  xlim([0 20])
  ylim([-3 3])
  grid on
  
  drawnow;
  frame = getframe(fig_anim);
  writeVideo(v,frame);
end
close(v);

fprintf('EXP dam DONE\n');

%% Local init function
function [h,u] = dam_init(x,Bc)
  % eta stratified:
  eta = zeros(size(x));
  eta(x<8) = 2;
  eta((x>=8) & (x<12)) = 1;
  eta(x>=12) = 1; % over step B=0.5
  
  h = max(0, eta - Bc);
  u = zeros(size(h));
end

function u = safe_u(H, HU)
  tol = 1e-10;
  u = zeros(size(H));
  wet = H > tol;
  u(wet) = HU(wet)./H(wet);
end
