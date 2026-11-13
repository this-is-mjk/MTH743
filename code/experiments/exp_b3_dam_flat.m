% exp_b3_dam_flat: basic idea 3 - the shallow water equations WITHOUT a bed: dam break against the exact Riemann solution
%   h_t + (hu)_x = 0, (hu)_t + (hu^2 + g h^2/2)_x = 0 on a flat bed.  Wet dam break (h_L = 2, h_R = 1) and dry-bed dam
%   break (h_L = 1, h_R = 0, Ritter's solution).  Domain [0,10], dam at x = 5, g = 9.81.
%   Schemes: Rusanov (local Lax-Friedrichs) first order, HLL first order (library, scheme 'hr' order 1 on a flat bed),
%   WENO-3 + SSP-RK3 (library, order 3).  Reference: LeVeque ch. 13 (exact solution), ch. 15 (approximate solvers).
set(0,'DefaultFigureVisible','off'); s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light';
root   = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end
addpath(fullfile(fileparts(mfilename('fullpath')),'..','lib'));

g = 9.81; x0 = 5; L = 10;
cases = struct('name',{'wet','dry'},'hL',{2,1},'hR',{1,0},'T',{0.8,0.6});
names = {'Rusanov (order 1)','HLL (order 1)','WENO-3 + SSP-RK3'};
col = [0.85 0.325 0.098; 0.466 0.674 0.188; 0 0.447 0.741];
fid = fopen(fullfile(resdir,'b3.txt'),'w');
Ns = [50 100 200 400 800];
for c = 1:2
  cs = cases(c);
  %% exact solution: structure of the wet case
  if c == 1
    [hs,us] = star_state(cs.hL,0,cs.hR,0,g);
    SR = us*hs/(hs-cs.hR);                                    % right shock speed (mass conservation across the shock)
    fprintf(fid,'Wet dam break h_L = %g, h_R = %g: star state h* = %.4f, u* = %.4f, right shock speed = %.4f, left rarefaction head speed = %.4f\n\n',cs.hL,cs.hR,hs,us,SR,-sqrt(g*cs.hL));
  else
    fprintf(fid,'Dry dam break h_L = %g, h_R = 0 (Ritter): front speed 2 sqrt(g h_L) = %.4f, rarefaction head speed -sqrt(g h_L) = %.4f\n\n',cs.hL,2*sqrt(g*cs.hL),-sqrt(g*cs.hL));
  end
  %% convergence
  E = zeros(3,numel(Ns)); minh = zeros(3,numel(Ns));
  for j = 1:numel(Ns)
    N = Ns(j); [xc,hex,uex] = exact_avg(cs,N,L,x0,g);
    for k = 1:3
      [h,u] = safe_solve(k,cs,N,L,x0,g,0.5);        % NaN if the run blows up (WENO-3 on a dry bed)
      E(k,j) = sum(abs(h-hex))*L/N; minh(k,j) = min(h);
    end
  end
  fprintf(fid,'%s bed: L1 error of h at t = %.1f (observed order)\n',cs.name,cs.T);
  fprintf(fid,'%-22s',''); fprintf(fid,'%18d',Ns); fprintf(fid,'\n');
  for k = 1:3
    fprintf(fid,'%-22s',names{k});
    for j = 1:numel(Ns)
      if isnan(E(k,j)), fprintf(fid,'%18s','blow-up');
      elseif j == 1 || isnan(E(k,j-1)), fprintf(fid,'%11.3e       ',E(k,j)); else, fprintf(fid,'%11.3e (%.2f)',E(k,j),log2(E(k,j-1)/E(k,j))); end
    end
    fprintf(fid,'\n');
  end
  fprintf(fid,'smallest depth at t = %.1f over all N (positivity): Rusanov %.3e;  HLL %.3e;  WENO-3 %.3e (NaN = run blew up)\n\n',cs.T,min(minh(1,:)),min(minh(2,:)),min(minh(3,:)));
  Eall{c} = E;
  %% profiles N = 200
  N = 200; [xc,hex,uex] = exact_avg(cs,N,L,x0,g); xe = linspace(0,L,4000); [he,ue] = exact_point(cs,xe,x0,g);
  fig = figure('Color','w','Position',[50 50 1400 560]); tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
  H = cell(3,1); U = cell(3,1); for k = 1:3, [H{k},U{k}] = safe_solve(k,cs,N,L,x0,g,0.5); end
  nexttile; hold on; box on; plot(xe,he,'Color',[0 0 0],'LineWidth',2,'DisplayName','exact');
  for k = 1:3, nm = names{k}; if all(isnan(H{k})), nm = [nm ' - blows up']; end
    plot(xc,H{k},'Color',col(k,:),'LineWidth',1.6,'Marker','o','MarkerSize',3,'LineStyle','-','DisplayName',nm); end
  set(gca,'FontSize',12); xlabel('x [m]'); ylabel('h [m]'); legend('Location','northeast'); xlim([0 L]); title('depth');
  nexttile; hold on; box on; plot(xe,ue,'Color',[0 0 0],'LineWidth',2,'DisplayName','exact');
  for k = 1:3, plot(xc,U{k},'Color',col(k,:),'LineWidth',1.6,'Marker','o','MarkerSize',3,'LineStyle','-','DisplayName',names{k}); end
  set(gca,'FontSize',12); xlabel('x [m]'); ylabel('u [m/s]'); xlim([0 L]); title('velocity');
  title(tl,sprintf('Dam break on a flat %s bed, h_L = %g, h_R = %g, t = %.1f s, N = %d',cs.name,cs.hL,cs.hR,cs.T,N),'FontSize',13);
  exportgraphics(fig,fullfile(imgdir,sprintf('b3_dam_%s.png',cs.name)),'Resolution',150); close(fig);
end
%% convergence figure
fig = figure('Color','w','Position',[50 50 1300 560]); tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); dxs = L./Ns;
for c = 1:2
  nexttile; hold on; box on; ls = {'-','--',':'};
  for k = 1:3, loglog(dxs,Eall{c}(k,:),'Color',col(k,:),'LineStyle',ls{k},'LineWidth',1.8,'Marker','o','MarkerSize',7,'DisplayName',names{k}); end
  loglog(dxs,Eall{c}(2,1)*(dxs/dxs(1)),'Color',[0 0 0],'LineStyle',':','LineWidth',1.2,'DisplayName','slope 1');
  set(gca,'XScale','log','YScale','log','FontSize',12); grid on; xlabel('\Delta x [m]'); ylabel('L^1 error of h'); legend('Location','southeast'); title([cases(c).name ' bed']);
end
title(tl,'Dam break on a flat bed: shocks limit all schemes to about first order in L^1','FontSize',13);
exportgraphics(fig,fullfile(imgdir,'b3_dam_convergence.png'),'Resolution',150); close(fig);
%% time step: first-order HLL scheme versus CFL number (wet case, N = 200)
cs = cases(1); N = 200; fprintf(fid,'Time-step limit of the first-order HLL scheme, wet dam break, N = %d, t = %.1f:\n',N,cs.T);
fprintf(fid,'%6s %10s %12s %12s\n','CFL','finished','max h','min h');
for cfl = [0.5 0.9 1.0 1.1 1.5]
  ok = true; try, [h,~] = solve(2,cs,N,L,x0,g,cfl); catch, ok = false; end
  if ok, fprintf(fid,'%6.2f %10s %12.4f %12.4f\n',cfl,'yes',max(h),min(h)); else, fprintf(fid,'%6.2f %10s %12s %12s\n',cfl,'blow-up','-','-'); end
end
fclose(fid); fprintf('EXP b3 DONE\n');

% ------------------------------------------------------------------------------------------------------------
function [h,u] = safe_solve(kind,cs,N,L,x0,g,cfl)
  try, [h,u] = solve(kind,cs,N,L,x0,g,cfl); catch, h = nan(1,N); u = nan(1,N); end
end
function [h,u] = solve(kind,cs,N,L,x0,g,cfl)
  init = @(x,Bc) deal(cs.hL*(x<x0) + cs.hR*(x>=x0), 0*x);
  if kind == 1, [h,u] = rusanov(cs,N,L,x0,g,cfl); return; end
  cfg = struct('N',N,'xl',0,'xr',L,'bed',make_bed('flat',[]),'init',init,'tend',cs.T,'tsave',cs.T,'bc','open','g',g,'scheme','hr','cfl',cfl);
  if kind == 2, cfg.order = 1; cfg.rk = 'euler'; else, cfg.order = 3; cfg.rk = 'ssprk3'; end
  out = swe1d_run(cfg); h = out.H{end}; hu = out.HU{end}; u = zeros(size(h)); w = h > 1e-10; u(w) = hu(w)./h(w);
end
function [h,u] = rusanov(cs,N,L,x0,g,cfl)
  dx = L/N; xc = ((1:N)-0.5)*dx; h = cs.hL*(xc<x0) + cs.hR*(xc>=x0); q = 0*h; t = 0;
  while t < cs.T - 1e-14
    w = h > 1e-10; u = zeros(size(h)); u(w) = q(w)./h(w); a = max(abs(u)+sqrt(g*h)); dt = min(cfl*dx/a,cs.T-t);
    hp = [h(1) h h(end)]; qp = [q(1) q q(end)]; up = zeros(size(hp)); wp = hp > 1e-10; up(wp) = qp(wp)./hp(wp);
    f1 = qp; f2 = qp.*up + 0.5*g*hp.^2; s = abs(up) + sqrt(g*hp); sm = max(s(1:end-1),s(2:end));
    F1 = 0.5*(f1(1:end-1)+f1(2:end)) - 0.5*sm.*(hp(2:end)-hp(1:end-1));
    F2 = 0.5*(f2(1:end-1)+f2(2:end)) - 0.5*sm.*(qp(2:end)-qp(1:end-1));
    h = h - dt/dx*(F1(2:end)-F1(1:end-1)); q = q - dt/dx*(F2(2:end)-F2(1:end-1)); t = t + dt;
  end
  u = zeros(size(h)); w = h > 1e-10; u(w) = q(w)./h(w);
end
function [xc,h,u] = exact_avg(cs,N,L,x0,g)       % exact cell averages (10 sample points per cell)
  nsub = 10; xs = ((0:N*nsub-1)+0.5)*(L/(N*nsub)); [hp,up] = exact_point(cs,xs,x0,g);
  h = mean(reshape(hp,nsub,N),1); u = mean(reshape(up,nsub,N),1); xc = ((1:N)-0.5)*(L/N);
end
function [hs,us] = star_state(hL,uL,hR,uR,g)
  fk = @(h,hk) (h>hk).*((h-hk).*sqrt(0.5*g*(1./h+1./hk))) + (h<=hk).*(2*(sqrt(g*h)-sqrt(g*hk)));
  F = @(h) fk(h,hL) + fk(h,hR) + (uR-uL);
  hs = fzero(F,[1e-8 10*max(hL,hR)]); us = 0.5*(uL+uR) + 0.5*(fk(hs,hR)-fk(hs,hL));
end
function [h,u] = exact_point(cs,x,x0,g)
  hL = cs.hL; hR = cs.hR; t = cs.T; S = (x-x0)/t; h = zeros(size(x)); u = zeros(size(x)); cL = sqrt(g*hL);
  if hR <= 0                                           % dry bed (Ritter)
    for i = 1:numel(x)
      if S(i) <= -cL, h(i) = hL; u(i) = 0;
      elseif S(i) <= 2*cL, c = (2*cL - S(i))/3; u(i) = S(i)+c; h(i) = c^2/g;
      else, h(i) = 0; u(i) = 0; end
    end
    return
  end
  [hs,us] = star_state(hL,0,hR,0,g); cstar = sqrt(g*hs);
  % left wave: rarefaction when hs < hL, shock when hs > hL;  right wave: shock (hs > hR) in the cases used here
  for i = 1:numel(x)
    Si = S(i);
    if hs <= hL                                        % left rarefaction
      if Si <= -cL, h(i) = hL; u(i) = 0;
      elseif Si <= us - cstar, c = (2*cL - Si)/3; u(i) = Si + c; h(i) = c^2/g;
      elseif Si <= us, h(i) = hs; u(i) = us;
      else, [h(i),u(i)] = right_side(Si,hs,us,hR); end
    else                                               % left shock
      SLs = -sqrt(g*hs*(hs+hL)/(2*hL));
      if Si <= SLs, h(i) = hL; u(i) = 0; elseif Si <= us, h(i) = hs; u(i) = us; else, [h(i),u(i)] = right_side(Si,hs,us,hR); end
    end
  end
end
function [h,u] = right_side(S,hs,us,hR)                % right of the contact (u = u*): shock, hs > hR
  SR = us*hs/(hs-hR);
  if S <= SR, h = hs; u = us; else, h = hR; u = 0; end
end
