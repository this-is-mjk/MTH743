% exp_b1_advection: basic idea 1 - what a scheme does to a transported profile (diffusion, dispersion, limiting)
%   u_t + a u_x = 0 on [0,1], periodic, a = 1, one period T = 1 (the exact solution is the initial profile).
%   Schemes: first-order upwind, Lax-Wendroff, MUSCL-minmod (SSP-RK2), WENO-3 (SSP-RK3), WENO-3 with forward Euler.
%   Initial data are exact CELL AVERAGES (so the orders of the high-order schemes are not polluted by the
%   point-value initialisation).  Reference: LeVeque ch. 4, 6, 8, 10; Mishra ch. 2, 5, 6.
set(0,'DefaultFigureVisible','off'); s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light';
root   = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end
addpath(fullfile(fileparts(mfilename('fullpath')),'..','lib'));

a = 1; T = 1;
smooth = @(x) exp(-((x-0.5)/0.1).^2);                                   % smooth pulse
mixed  = @(x) exp(-((x-0.3)/0.07).^2) + 1.0*(x>0.6 & x<0.8);           % Gaussian + square wave
schemes = {'upwind','laxwendroff','muscl','weno3','weno3euler'};
names   = {'upwind (order 1)','Lax-Wendroff (order 2)','MUSCL-minmod + SSP-RK2','WENO-3 + SSP-RK3','WENO-3 + forward Euler'};
cfls    = [0.8 0.8 0.4 0.4 0.4];
col = [0.85 0.325 0.098; 0.494 0.184 0.556; 0.466 0.674 0.188; 0 0.447 0.741; 0.3 0.3 0.3];

fid = fopen(fullfile(resdir,'b1.txt'),'w');
%% 1. convergence on the smooth pulse
Ns = [40 80 160 320 640]; E = zeros(numel(schemes),numel(Ns));
for k = 1:numel(schemes)
  for j = 1:numel(Ns)
    u0 = cellavg(smooth,Ns(j)); u = advect(u0,Ns(j),a,T,schemes{k},cfls(k));
    E(k,j) = sum(abs(u-u0))/Ns(j);
  end
end
fprintf(fid,'Smooth pulse exp(-((x-0.5)/0.1)^2), one period: L1 error (observed order)\n');
fprintf(fid,'%-28s','scheme'); fprintf(fid,'%16d',Ns); fprintf(fid,'\n');
for k = 1:numel(schemes)
  fprintf(fid,'%-28s',names{k});
  for j = 1:numel(Ns)
    if j == 1 || ~isfinite(E(k,j)) || ~isfinite(E(k,j-1)) || E(k,j) == 0, fprintf(fid,'%11.3e      ',E(k,j));
    else, fprintf(fid,'%11.3e (%.2f)',E(k,j),log2(E(k,j-1)/E(k,j))); end
  end
  fprintf(fid,'\n');
end
%% 2. mixed profile: shape of the solution after 1 and 5 periods (N = 200)
N = 200; x = ((1:N)-0.5)/N; fig = figure('Color','w','Position',[50 50 1300 800]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
fprintf(fid,'\nGaussian + square wave, N = %d: L1 error, min, max, total variation (initial TV = %.3f, range [0,1.0])\n',N,tv(cellavg(mixed,N)));
fprintf(fid,'%-28s %-8s %12s %10s %10s %10s\n','scheme','periods','L1 error','min','max','TV');
u0 = cellavg(mixed,N);
for k = 1:numel(schemes)
  for P = [1 5]
    u = advect(u0,N,a,P*T,schemes{k},cfls(k));
    fprintf(fid,'%-28s %-8d %12.3e %10.3f %10.3f %10.3f\n',names{k},P,sum(abs(u-u0))/N,min(u),max(u),tv(u));
    if P == 1 && k <= 4
      nexttile; hold on; box on;
      plot(x,u0,'Color',[0 0 0],'LineStyle','-','LineWidth',1.6);
      plot(x,u,'Color',col(k,:),'LineStyle','-','LineWidth',1.8,'Marker','o','MarkerSize',3);
      xlim([0 1]); ylim([-0.35 1.4]); set(gca,'FontSize',12); xlabel('x'); ylabel('u');
      title(sprintf('%s: L^1 = %.2e, min = %.2f, max = %.2f',names{k},sum(abs(u-u0))/N,min(u),max(u)),'FontSize',11);
      legend({'exact','numerical'},'Location','northeast');
    end
  end
end
title(tl,'One period of linear advection, N = 200 (exact = initial profile)','FontSize',13);
exportgraphics(fig,fullfile(imgdir,'b1_advection_profiles.png'),'Resolution',150); close(fig);
%% 3. convergence figure
fig = figure('Color','w','Position',[50 50 900 600]); hold on; box on; dxs = 1./Ns; ls = {'-','--','-.',':','-'};
for k = 1:4
  loglog(dxs,E(k,:),'Color',col(k,:),'LineStyle',ls{k},'LineWidth',1.8,'Marker','o','MarkerSize',7,'DisplayName',names{k});
end
for p = 1:3
  loglog(dxs,E(1,1)*(dxs/dxs(1)).^p*0.25,'Color',[0 0 0],'LineStyle',':','LineWidth',1.2,'HandleVisibility','off');
  text(dxs(end)*1.15,E(1,1)*0.25*(dxs(end)/dxs(1))^p,sprintf('slope %d',p),'FontSize',11);
end
set(gca,'XScale','log','YScale','log','FontSize',12); grid on; xlabel('\Delta x'); ylabel('L^1 error after one period');
title('Smooth pulse: order of the schemes'); legend('Location','southeast'); xlim([dxs(end)*0.8 dxs(1)*1.5]);
exportgraphics(fig,fullfile(imgdir,'b1_advection_convergence.png'),'Resolution',150); close(fig);
fclose(fid); fprintf('EXP b1 DONE\n');

% ------------------------------------------------------------------------------------------------------------
function u0 = cellavg(f,N)
  nsub = 200; dx = 1/N; xs = ((0:N*nsub-1)+0.5)/(N*nsub);
  u0 = mean(reshape(f(xs),nsub,N),1);
end
function t = tv(u), t = sum(abs(diff([u u(1)]))); end
function u = advect(u,N,a,T,scheme,cfl)
  dx = 1/N; nst = ceil(T/(cfl*dx/a)); dt = T/nst; lam = a*dt/dx;
  for n = 1:nst
    switch scheme
      case 'upwind',      u = u - lam*(u - circshift(u,1));
      case 'laxwendroff', up = circshift(u,-1); um = circshift(u,1); u = u - lam/2*(up-um) + lam^2/2*(up-2*u+um);
      case 'muscl',       u1 = u + dt*rhs(u,a,dx,'muscl'); u = 0.5*u + 0.5*(u1 + dt*rhs(u1,a,dx,'muscl'));
      case 'weno3',       u1 = u + dt*rhs(u,a,dx,'weno3'); u2 = 0.75*u + 0.25*(u1 + dt*rhs(u1,a,dx,'weno3'));
                          u = u/3 + 2/3*(u2 + dt*rhs(u2,a,dx,'weno3'));
      case 'weno3euler',  u = u + dt*rhs(u,a,dx,'weno3');
    end
    if ~all(isfinite(u)) || max(abs(u)) > 1e6, u = nan(size(u)); return; end
  end
end
function r = rhs(u,a,dx,kind)
  N = numel(u);
  switch kind
    case 'muscl'
      d1 = u - circshift(u,1); d2 = circshift(u,-1) - u;
      sl = (sign(d1)+sign(d2))/2.*min(abs(d1),abs(d2));          % minmod
      uf = u + sl/2;                                             % value at the right face of each cell
      F = a*uf; r = -(F - circshift(F,1))/dx;
    case 'weno3'
      up = [u(end-1:end) u u(1:2)]; vL = weno3_faces(up); F = a*vL;
      r = -(F(2:N+1) - F(1:N))/dx;
  end
end
