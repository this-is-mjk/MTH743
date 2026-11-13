% exp_b4_balance_law: basic idea 4 - what "well balanced" means, in the simplest balance law
%   u_t + u_x = -beta u  on [0,1], inflow u(0,t) = 1, beta = 1.   Steady state: u(x) = exp(-beta x)
%   (transport and source cancel exactly, like pressure and bed slope in a lake at rest).
%   NAIVE scheme (transport by upwind, source at the cell):  u_i <- u_i - lam (u_i - u_{i-1}) - dt beta u_i
%       its own steady state is u_i = u_{i-1}/(1 + beta dx), not exp(-beta dx) u_{i-1}  -> error O(dx).
%   BALANCED scheme (source integrated along the transport):   u_i <- u_i - lam (u_i - exp(-beta dx) u_{i-1})
%       steady state u_i = exp(-beta dx) u_{i-1}: exact.  Both are first order accurate and consistent.
%   Reference: LeVeque sec. 17.8 / 17.14 (steady states, quasi-steady problems).
set(0,'DefaultFigureVisible','off'); s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light';
root   = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

beta = 1; lam = 0.8; T = 5;          % T = 5 transit times: the initial data have been replaced several times
ustar = @(x) exp(-beta*x);
col = [0.85 0.325 0.098; 0 0.447 0.741];
fid = fopen(fullfile(resdir,'b4.txt'),'w');
%% 1. start from the exact steady state: does the scheme keep it?
Ns = [25 50 100 200 400 800];
fprintf(fid,'Steady state u = exp(-x), start from the exact steady state, run to t = %g\n',T);
fprintf(fid,'%6s %18s %18s %10s %10s\n','N','naive max error','balanced max error','naive ord.','ratio');
Enaive = zeros(size(Ns)); Ebal = zeros(size(Ns));
for j = 1:numel(Ns)
  N = Ns(j); xc = ((1:N)-0.5)/N; u0 = cellavg_exp(N,beta);
  Enaive(j) = max(abs(run(u0,N,beta,lam,T,'naive') - u0)); Ebal(j) = max(abs(run(u0,N,beta,lam,T,'balanced') - u0));
  if j == 1, o = NaN; else, o = log2(Enaive(j-1)/Enaive(j)); end
  fprintf(fid,'%6d %18.3e %18.3e %10.2f %10.1e\n',N,Enaive(j),Ebal(j),o,Enaive(j)/max(Ebal(j),eps));
end
%% 2. small perturbation travelling on the steady state (the tsunami situation)
N = 100; dx = 1/N; xc = ((1:N)-0.5)/N; eps0 = 1e-3; tp = 0.4;
pert = @(x) eps0*exp(-((x-0.2)/0.05).^2);
u0 = cellavg_exp(N,beta) + cellavg_fun(pert,N);
exact = @(x) ustar(x) + eps0*exp(-beta*tp)*exp(-((x-tp-0.2)/0.05).^2);   % linear equation: the pulse is transported and damped
uex = cellavg_fun(exact,N);
un = run(u0,N,beta,lam,tp,'naive'); ub = run(u0,N,beta,lam,tp,'balanced');
sigE = @(u) u - cellavg_exp(N,beta);                       % deviation from the steady state = signal + scheme error
fprintf(fid,'\nPerturbation of size %.0e travelling on the steady state, N = %d, t = %.1f (signal amplitude at t: %.2e)\n',eps0,N,tp,eps0*exp(-beta*tp));
fprintf(fid,'  naive:    max |u - exact| = %.3e   (relative to the signal amplitude: %.1f)\n',max(abs(un-uex)),max(abs(un-uex))/(eps0*exp(-beta*tp)));
fprintf(fid,'  balanced: max |u - exact| = %.3e   (relative to the signal amplitude: %.2f)\n',max(abs(ub-uex)),max(abs(ub-uex))/(eps0*exp(-beta*tp)));
fclose(fid);
%% figures
fig = figure('Color','w','Position',[50 50 1500 520]); tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on; N1 = 50; xc1 = ((1:N1)-0.5)/N1; u01 = cellavg_exp(N1,beta);
plot(xc1,run(u01,N1,beta,lam,T,'naive')-u01,'Color',col(1,:),'LineWidth',1.8,'Marker','o','MarkerSize',4,'DisplayName','naive');
plot(xc1,run(u01,N1,beta,lam,T,'balanced')-u01,'Color',col(2,:),'LineWidth',1.8,'Marker','s','MarkerSize',4,'DisplayName','balanced');
set(gca,'FontSize',12); xlabel('x'); ylabel('u - u_{steady}'); legend('Location','southwest'); title('error after t = 5, N = 50','FontSize',11);
nexttile; hold on; box on; dxs = 1./Ns;
loglog(dxs,Enaive,'Color',col(1,:),'LineWidth',1.8,'Marker','o','MarkerSize',7,'DisplayName','naive');
loglog(dxs,max(Ebal,1e-17),'Color',col(2,:),'LineWidth',1.8,'Marker','s','MarkerSize',7,'LineStyle','--','DisplayName','balanced');
loglog(dxs,Enaive(1)*dxs/dxs(1),'Color',[0 0 0],'LineStyle',':','LineWidth',1.2,'DisplayName','slope 1');
set(gca,'XScale','log','YScale','log','FontSize',12); grid on; xlabel('\Delta x'); ylabel('max error'); legend('Location','east'); title('steady-state error versus \Delta x'); ylim([1e-17 1]);
nexttile; hold on; box on;
plot(xc,uex-cellavg_exp(N,beta),'Color',[0 0 0],'LineWidth',2,'DisplayName','exact signal');
plot(xc,un-cellavg_exp(N,beta),'Color',col(1,:),'LineWidth',1.8,'Marker','o','MarkerSize',3,'DisplayName','naive');
plot(xc,ub-cellavg_exp(N,beta),'Color',col(2,:),'LineWidth',1.8,'Marker','s','MarkerSize',3,'DisplayName','balanced');
set(gca,'FontSize',12); xlabel('x'); ylabel('u - u_{steady}'); legend('Location','northeast'); title('perturbation 10^{-3} at t = 0.4, N = 100');
title(tl,'Scalar balance law u_t + u_x = -u: the naive scheme drifts away from the steady state and hides small waves','FontSize',13);
exportgraphics(fig,fullfile(imgdir,'b4_balance_law.png'),'Resolution',150); close(fig);
fprintf('EXP b4 DONE\n');

% ------------------------------------------------------------------------------------------------------------
function u = run(u,N,beta,lam,T,kind)
  dx = 1/N; nst = ceil(T/(lam*dx)); dt = T/nst; l = dt/dx; ex = exp(-beta*dx);
  for n = 1:nst
    um = [1 u(1:end-1)];                             % inflow boundary value u(0) = 1 (cell 0 holds the exact steady inflow)
    um(1) = (exp(beta*dx) - 1)/(beta*dx);            % exact average of exp(-beta x) over the ghost cell [-dx,0]: inflow = the steady state
    switch kind
      case 'naive',    u = u - l*(u - um) - dt*beta*u;
      case 'balanced', u = u - l*(u - ex*um);
    end
  end
end
function u0 = cellavg_exp(N,beta)
  dx = 1/N; xe = (0:N)*dx; u0 = (exp(-beta*xe(1:end-1)) - exp(-beta*xe(2:end)))/(beta*dx);   % exact cell averages of exp(-beta x)
end
function u0 = cellavg_fun(f,N)
  nsub = 100; xs = ((0:N*nsub-1)+0.5)/(N*nsub); u0 = mean(reshape(f(xs),nsub,N),1);
end
