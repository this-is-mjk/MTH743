% exp_p2: Phase 2, hydrostatic reconstruction on a lake at rest over three beds; small-wave animation
function exp_p2
set(0,'DefaultFigureVisible','off');
s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light';
root = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

g = 9.81;
xl=0; xr=10; eta0=1;

% Beds definition
beds = { ...
  struct('name','gaussian','type','gaussian','params',[0.4 5 0.5]), ...
  struct('name','cosbump' ,'type','cosbump' ,'params',[0.5 5 2]), ...
  struct('name','gaussian2','type','gaussian','params',[0.8 3 0.4]) ...
  };

Nvals = [50 100 200 400];
schemes = {'central','hr','hr'};
orders = [1 1 3];
rkmap = {'euler','euler','ssprk3'};
cflmap = [0.5 0.5 0.4];
c = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556; 0.301 0.745 0.933; 0.2 0.2 0.2];

% Prepare results storage for printing
results = cell(0,5);
% Format: bed, N, scheme-order, max|u|, max|eta-eta0|

% LOOP over beds and N and schemes
for ib = 1:numel(beds)
  bedstruct = beds{ib};
  bed = make_bed(bedstruct.type, bedstruct.params);
  for iN = 1:numel(Nvals)
    N = Nvals(iN);
    for isch = 1:numel(schemes)
      scheme = schemes{isch};
      order = orders(isch);
      rk = rkmap{isch};
      cfl = cflmap(isch);

      % init lake at rest
      init = @(x,Bc) deal( max(0, eta0 - Bc), zeros(size(x)) );
      
      cfg = struct();
      cfg.N = N;
      cfg.xl = xl;
      cfg.xr = xr;
      cfg.bed = bed;
      cfg.init = init;
      cfg.tend = 10;
      cfg.scheme = scheme;
      cfg.order = order;
      cfg.rk = rk;
      cfg.bc = 'wall';
      cfg.cfl = cfl;
      cfg.monitor = true;
      % Save outputs only at final time
      cfg.tsave = 10;
      
      out = swe1d_run(cfg);
      
      % Compute max|u| over whole run from monitor data
      maxu = max(out.mon.umax);
      
      % Compute max|eta-eta0| at final time:
      % eta = h + B cell center
      h_final = out.H{1};
      Bc = get_bc_bed(bed, N, xl, xr, order, scheme, 'wall'); % helper to get bed at centers
      eta_final = h_final + Bc;
      max_eta_diff = max(abs(eta_final - eta0));
      
      results(end+1,:) = {bedstruct.name, N, sprintf('%s-%d',scheme,order), maxu, max_eta_diff};

    end
  end
end

% Save results text file and print
fname = fullfile(resdir,'p2.txt');
fid = fopen(fname,'w');
fprintf(fid, 'bed\tN\tscheme-order\tmax|u|\tmax|eta-eta0|\n');
fprintf('bed\tN\tscheme-order\tmax|u|\tmax|eta-eta0|\n');
for k=1:size(results,1)
  r = results(k,:);
  fprintf(fid,'%s\t%d\t%s\t%g\t%g\n', r{1}, r{2}, r{3}, r{4}, r{5});
  fprintf('%s\t%d\t%s\t%g\t%g\n', r{1}, r{2}, r{3}, r{4}, r{5});
end

% === FIGURE 1: semilogy max|u| vs N for all beds and schemes, 3 panels ===
fig1 = figure('Color','w','Position',[100 100 900 600]);
for ib=1:numel(beds)
  subplot(1,3,ib);
  hold on
  for isch = 1:numel(schemes)
    sname = sprintf('%s-%d',schemes{isch},orders(isch));
    maxuN = nan(1,numel(Nvals));
    for iN=1:numel(Nvals)
      mask = strcmp(results(:,1),beds{ib}.name) & ...
             [results{:,2}]'==Nvals(iN) & ...
             strcmp(results(:,3),sname);
      if any(mask)
        maxuN(iN) = results{find(mask,1),4};
      end
    end
    plot(Nvals, maxuN, '-o', 'DisplayName', sname, 'LineWidth',1.5, 'Color', c(isch,:));
  end
  % machine epsilon line
  yline(eps,'k--','LineWidth',1.5,'DisplayName','machine eps');
  set(gca,'YScale','log');
  xlabel('N (cells)');
  ylabel('max |u|');
  title(sprintf('Bed: %s',beds{ib}.name));
  legend('Location','best','FontSize',10);
  grid on; box on;
end
set(findall(fig1,'-property','FontSize'),'FontSize',11);
exportgraphics(fig1, fullfile(imgdir,'p2_maxu_vs_N.png'), 'Resolution',150);

% === FIGURE 2: for N=200, bed (i) gaussian, plot u(x) at t=10 for central-1 and hr-1 ===
N = 200; ib = 1;
bedstruct = beds{ib}; 
bed = make_bed(bedstruct.type, bedstruct.params);
tend = 10;

% Run central-1
cfg1 = struct('N',N,'xl',xl,'xr',xr,'bed',bed,'init',@(x,Bc) deal(max(0,eta0-Bc),zeros(size(x))),...
  'tend',tend,'scheme','central','order',1,'rk','euler','bc','wall','cfl',0.5,'tsave',tend);
out1 = swe1d_run(cfg1);
u1 = out1.HU{1} ./ out1.H{1};
maxu1 = max(abs(u1));

% Run hr-1
cfg2 = struct('N',N,'xl',xl,'xr',xr,'bed',bed,'init',@(x,Bc) deal(max(0,eta0-Bc),zeros(size(x))),...
  'tend',tend,'scheme','hr','order',1,'rk','euler','bc','wall','cfl',0.5,'tsave',tend);
out2 = swe1d_run(cfg2);
u2 = out2.HU{1} ./ out2.H{1};
maxu2 = max(abs(u2));

xf = linspace(xl,xr,N); % cell centers approx as out.x

fig2 = figure('Color','w','Position',[100 100 900 600]);
subplot(1,2,1)
plot(out1.x,u1,'Color', c(1,:),'LineWidth',1.5);
ylim([-1 1]*maxu1*1.1);
xlabel('x (m)');
ylabel('u (m/s)');
title(sprintf('central-1, max|u|=%.2e',maxu1));
grid on; box on;

subplot(1,2,2)
plot(out2.x,u2,'Color', c(2,:),'LineWidth',1.5);
ylim([-1 1]*maxu2*1.1);
xlabel('x (m)');
ylabel('u (m/s)');
title(sprintf('hr-1, max|u|=%.2e',maxu2));
grid on; box on;

set(findall(fig2,'-property','FontSize'),'FontSize',11);
exportgraphics(fig2, fullfile(imgdir,'p2_u_compare_N200_bed1.png'), 'Resolution',150);

% === Perturbation robustness: N=200, bed(i), hr-3, perturb eta0 by small Gaussian ===
N = 200; ib = 1;
bedstruct = beds{ib};
bed = make_bed(bedstruct.type, bedstruct.params);
tend = 2;
xg = linspace(xl,xr,N);

% Initial condition eta=1+1e-3*exp(-(x-2).^2/0.05), u=0
init_pert = @(x,Bc) deal(max(0,eta0 + 1e-3*exp(-(x-2).^2/0.05) - Bc), zeros(size(x)));

cfgp = struct('N',N,'xl',xl,'xr',xr,'bed',bed,'init',init_pert,'tend',tend,'scheme','hr','order',3,...
  'rk','ssprk3','bc','wall','cfl',0.4,'tsave',linspace(0,tend,101));
outp = swe1d_run(cfgp);

% Compute eta-1 at chosen times:
eta_t = cellfun(@(h) h + get_bc_bed(bed,N,xl,xr,3,'hr','wall'), outp.H, 'UniformOutput', false);
eta1 = cellfun(@(eta) eta - 1, eta_t, 'UniformOutput', false);

% Plot at t=0,0.5,1,2
tplot = [0 0.5 1 2];
[~,indts] = arrayfun(@(t) min(abs(outp.t - t)), tplot);

fig3 = figure('Color','w','Position',[100 100 900 600]);
for i=1:numel(tplot)
  subplot(2,2,i);
  plot(outp.x,eta1{indts(i)}, 'Color', c(1,:),'LineWidth',1.5);
  xlabel('x (m)');
  ylabel('\eta - 1 (m)');
  title(sprintf('t=%.1f s', tplot(i)));
  grid on; box on;
end
sgtitle('Perturbation \eta-1 vs x at selected times (bed gaussian, hr-3, N=200)');

set(findall(fig3,'-property','FontSize'),'FontSize',11);
exportgraphics(fig3, fullfile(imgdir,'p2_perturbation_eta.png'), 'Resolution',150);

% Max |eta-1| away from pulse region x \in [0,10] excluding [1.75,2.25] where pulse is
ex_area_mask = (outp.x < 1.75) | (outp.x > 2.25);
max_eta_away = zeros(size(outp.t));
for i=1:numel(outp.t)
  max_eta_away(i) = max(abs(eta1{i}(ex_area_mask)));
end

% Print max|eta-1| away from pulse at final time t=2
fid = fopen(fname,'a');
fprintf(fid, '\nPerturbation robustness run (N=200, bed %s, hr-3):\n', bedstruct.name);
fprintf(fid, 'max|eta-1| away from pulse at times t=0,0.5,1,2: %g %g %g %g\n', max_eta_away(indts));
fprintf('Perturbation robustness max|eta-1| away from pulse: %g at t=2\n', max_eta_away(indts(end)));
fclose(fid);

% === Create MP4 of perturbation run ===
v = VideoWriter(fullfile(imgdir,'p2_perturbation.mp4'),'MPEG-4');
v.FrameRate = 20;
open(v);
fig4 = figure('Color','w','Position',[100 100 900 600]);
for k=1:numel(outp.t)
  plot(outp.x, eta1{k}, 'Color', c(1,:),'LineWidth',1.5);
  ylim([-2e-3 2e-3]);
  xlabel('x (m)'); ylabel('\eta - 1 (m)');
  title(sprintf('Perturbation at t = %.3f s', outp.t(k)));
  grid on; box on;
  drawnow;
  frame = getframe(fig4);
  writeVideo(v,frame);
end
close(v);
close(fig4);

% === FINAL line worst max|u| of each method ===
% Gather max|u| per method (scheme-order)
mnames = unique(results(:,3));
worstmaxu = zeros(size(mnames));
for i=1:numel(mnames)
  sel = strcmp(results(:,3),mnames{i});
  worstmaxu(i) = max(cell2mat(results(sel,4)));
end
fid = fopen(fname,'a');
fprintf(fid, '\nWorst max|u| per method:\n');
for i=1:numel(mnames)
  fprintf(fid, '%s: %.3e\n', mnames{i}, worstmaxu(i));
  fprintf('%s: %.3e\n', mnames{i}, worstmaxu(i));
end
fclose(fid);

fprintf('EXP p2 DONE\n');

end

% helper function to get bed at cell centers with proper padding
function Bc = get_bc_bed(bed, N, xl, xr, order, scheme, bc)
  % Reproduce how swe1d_setup computes P.Bc with padding at order and scheme.
  cfg.N = N; cfg.xl=xl; cfg.xr=xr; cfg.bed=bed; cfg.scheme=scheme; cfg.order=order; cfg.bc=bc;
  P = swe1d_setup(cfg);
  % Return only cell central bed (physical cells only)
  Bc = P.Bc(3:N+2);
end
