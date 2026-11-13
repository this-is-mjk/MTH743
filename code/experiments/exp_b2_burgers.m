% exp_b2_burgers: basic idea 2 - conservation form decides the speed of a shock
%   Burgers u_t + (u^2/2)_x = 0 with the Riemann data u_L = 1, u_R = 0 (shock, Rankine-Hugoniot speed (u_L+u_R)/2 = 1/2)
%   and u_L = 0, u_R = 1 (rarefaction fan u = x'/t).   Domain [0,1], jump at x = 0.3, t = 0.4  =>  exact shock at x = 0.5.
%   Schemes: Godunov (conservative), Rusanov (conservative), and the NON-conservative upwind scheme for u_t + u u_x = 0.
%   Reference: LeVeque sec. 12.9 ("the importance of conservation form") and the Lax-Wendroff theorem (sec. 12.10).
set(0,'DefaultFigureVisible','off'); s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light';
root   = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

x0 = 0.3; T = 0.4; cfl = 0.8;
fid = fopen(fullfile(resdir,'b2.txt'),'w');
names = {'Godunov (conservative)','Rusanov (conservative)','upwind u_t + u u_x = 0 (non-conservative)'};
kinds = {'godunov','rusanov','nonconservative'};
col = [0 0.447 0.741; 0.466 0.674 0.188; 0.85 0.325 0.098];
%% shock: position and mass
Ns = [50 100 200 400 800];
fprintf(fid,'Burgers shock u_L = 1, u_R = 0, exact shock position at t = %.1f: x = %.3f (speed 1/2)\n',T,x0+0.5*T);
fprintf(fid,'%-44s %6s %12s %12s %14s\n','scheme','N','shock x','x error','L1 error');
for k = 1:3
  for N = Ns
    [xc,u] = burgers(1,0,N,x0,T,cfl,kinds{k}); ex = double(xc < x0+0.5*T);
    xs = shockpos(xc,u);
    fprintf(fid,'%-44s %6d %12.4f %12.4f %14.4e\n',names{k},N,xs,xs-(x0+0.5*T),sum(abs(u-ex))/N);
  end
end
%% mass balance for the shock run N = 200 (inflow 1 enters at x = 0, nothing leaves at x = 1 before the shock arrives)
N = 200; fprintf(fid,'\nMass check, N = %d: integral of u over [0,1] at t = %.1f; exact = %.3f (initial 0.3 + inflow flux f(1)*t = 0.2)\n',N,T,x0+0.5*T);
for k = 1:3
  [xc,u] = burgers(1,0,N,x0,T,cfl,kinds{k}); fprintf(fid,'%-44s integral u = %.4f\n',names{k},sum(u)/N);
end
%% rarefaction: entropy-satisfying fan
N = 200; fig = figure('Color','w','Position',[50 50 1400 560]); tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on;
xe = linspace(0,1,2000); plot(xe,double(xe < x0+0.5*T),'Color',[0 0 0],'LineWidth',2,'DisplayName','exact shock');
for k = 1:3
  [xc,u] = burgers(1,0,N,x0,T,cfl,kinds{k}); plot(xc,u,'Color',col(k,:),'LineWidth',1.8,'Marker','o','MarkerSize',3,'LineStyle','-','DisplayName',names{k});
end
xlim([0 1]); ylim([-0.1 1.2]); set(gca,'FontSize',12); xlabel('x'); ylabel('u'); legend('Location','northeast'); title(sprintf('shock u_L=1, u_R=0 at t=%.1f, N=%d',T,N));
nexttile; hold on; box on;
fan = @(x) min(1,max(0,(x-x0)/T)); plot(xe,fan(xe),'Color',[0 0 0],'LineWidth',2,'DisplayName','exact fan');
for k = 1:3
  [xc,u] = burgers(0,1,N,x0,T,cfl,kinds{k}); plot(xc,u,'Color',col(k,:),'LineWidth',1.8,'Marker','o','MarkerSize',3,'LineStyle','-','DisplayName',names{k});
  fprintf(fid,'rarefaction u_L=0,u_R=1, N=%d, %-44s L1 error = %.3e\n',N,names{k},sum(abs(u-fan(xc)))/N);
end
xlim([0 1]); ylim([-0.1 1.2]); set(gca,'FontSize',12); xlabel('x'); ylabel('u'); legend('Location','northwest'); title(sprintf('rarefaction u_L=0, u_R=1 at t=%.1f, N=%d',T,N));
title(tl,'Burgers equation: the conservative schemes move the shock at the right speed, the non-conservative one does not','FontSize',13);
exportgraphics(fig,fullfile(imgdir,'b2_burgers.png'),'Resolution',150); close(fig);
fclose(fid); fprintf('EXP b2 DONE\n');

% ------------------------------------------------------------------------------------------------------------
function [xc,u] = burgers(uL,uR,N,x0,T,cfl,kind)
  dx = 1/N; xc = ((1:N)-0.5)*dx; u = uL*(xc<x0) + uR*(xc>=x0);
  t = 0;
  while t < T - 1e-14
    dt = min(cfl*dx/max(abs(u)+1e-12),T-t); up = [u(1) u u(end)];   % transmissive ghost cells
    switch kind
      case 'godunov',  F = godunov_flux(up(1:end-1),up(2:end)); u = u - dt/dx*(F(2:end)-F(1:end-1));
      case 'rusanov',  a = max(abs(up(1:end-1)),abs(up(2:end))); F = 0.5*(up(1:end-1).^2/2+up(2:end).^2/2) - 0.5*a.*(up(2:end)-up(1:end-1));
                       u = u - dt/dx*(F(2:end)-F(1:end-1));
      case 'nonconservative'    % upwind discretisation of u_t + u u_x = 0 (u >= 0): u_i - dt/dx u_i (u_i - u_{i-1})
                       u = u - dt/dx*u.*(u - up(1:end-2));
    end
    t = t + dt;
  end
end
function F = godunov_flux(ul,ur)       % exact Riemann flux for the convex flux f = u^2/2
  f = @(v) v.^2/2; F = zeros(size(ul));
  sh = ul > ur;  s = (ul+ur)/2;                       % shock: take the left state if s > 0
  F(sh) = f(ul(sh)).*(s(sh) > 0) + f(ur(sh)).*(s(sh) <= 0);
  rf = ~sh;                                           % rarefaction: minimum of f over [ul,ur], 0 if it contains 0
  F(rf) = min(f(ul(rf)),f(ur(rf))); F(rf & ul <= 0 & ur >= 0) = 0;
end
function xs = shockpos(xc,u)           % position where u crosses 1/2 (first crossing from the left)
  i = find(u(1:end-1) >= 0.5 & u(2:end) < 0.5,1,'first');
  if isempty(i), xs = NaN; else, xs = xc(i) + (u(i)-0.5)/(u(i)-u(i+1))*(xc(i+1)-xc(i)); end
end
