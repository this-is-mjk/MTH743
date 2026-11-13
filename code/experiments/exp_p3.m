% exp_p3: Phase 3, WENO-3 + SSP-RK3: accuracy, cost, narrow pulse, well-balancing
function exp_p3
set(0,'DefaultFigureVisible','off'); 
s=settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue='light';

root = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

% Constants, parameters
g = 9.81;

%% 1. Smooth accuracy test periodic BC (bc='periodic')
% bed = sin2, init h = 5+exp(cos(2pi x)), hu=sin(cos(2pi x))
bc = 'periodic'; 
bed = make_bed('sin2',[]);
tend = 0.1;

Nvec = [25 50 100 200 400];
% Reference: hr, order 3, ssprk3, N=3200, cfl=0.2
Nref = 3200; cfl_ref = 0.2;
cfgref.N = Nref; cfgref.xl=0; cfgref.xr=1; cfgref.bc=bc;
cfgref.bed = bed; cfgref.tend = tend;
cfgref.scheme = 'hr'; 
cfgref.order = 3;
cfgref.rk = 'ssprk3';
cfgref.cfl = cfl_ref;

initfun = @(x,Bc) deal(5 + exp(cos(2*pi*x)), sin(cos(2*pi*x))./(5 + exp(cos(2*pi*x))));
cfgref.init = initfun;
cfgref.tsave = tend;

outref = swe1d_run(cfgref);

% Exact cell averaging of ref solution for each coarse N
xcell_mids = cell(size(Nvec));
href_coarse = cell(size(Nvec));
huref_coarse = cell(size(Nvec));
dxref = 1/Nref;

for i=1:numel(Nvec)
    N = Nvec(i);
    dx = 1/N;
    % Cells indices in ref for each coarse cell
    bins = reshape(1:Nref, [], Nref/N);
    % mean h and hu using reshaping
    h_ref = outref.H{1}; hu_ref = outref.HU{1};
    href_coarse{i} = mean(reshape(h_ref, Nref/N, N),1);
    huref_coarse{i} = mean(reshape(hu_ref, Nref/N, N),1);
end

% Compute coarse solutions for (order1,euler,cfl0.5) and (order3,ssprk3,cfl0.2)
methods = { struct('order',1,'rk','euler','cfl',0.5), ...
            struct('order',3,'rk','ssprk3','cfl',0.2) };

L1h = zeros(numel(methods),numel(Nvec));
L1hu = zeros(numel(methods),numel(Nvec));

for m=1:numel(methods)
    method = methods{m};
    for i=1:numel(Nvec)
        N = Nvec(i);
        dx = 1/N;
        cfg.N = N;
        cfg.xl = 0; cfg.xr = 1;
        cfg.bc = bc;
        cfg.bed = bed;
        cfg.tend = tend;
        cfg.scheme = 'hr';
        cfg.order = method.order;
        cfg.rk = method.rk;
        cfg.cfl = method.cfl;
        cfg.init = initfun;
        cfg.tsave = tend;

        out = swe1d_run(cfg);

        hnum = out.H{1};
        hunum = out.HU{1};

        % Compute L1 error vs ref averaged on coarse grid
        href = href_coarse{i};
        huref = huref_coarse{i};

        L1h(m,i) = sum(abs(hnum - href))*dx;
        L1hu(m,i) = sum(abs(hunum - huref))*dx;
    end
end

% Compute observed order p for h and hu for each method, between successive N
ord_h = log2(L1h(:,1:end-1)./L1h(:,2:end));
ord_hu = log2(L1hu(:,1:end-1)./L1hu(:,2:end));

% Print + save table to text file
fname1 = fullfile(resdir,'p3_smooth_accuracy.txt');
fid = fopen(fname1,'w');

headers = {'N','L1(h) order1', 'order1 order', 'L1(hu) order1', 'order1 order', ...
           'L1(h) order3', 'order3 order', 'L1(hu) order3', 'order3 order'};
fprintf(fid, '%10s %14s %8s %14s %8s %14s %8s %14s %8s\n', headers{:});
fprintf( ...
    ' %10s %14s %8s %14s %8s %14s %8s %14s %8s\n', ...
    '----','--------------','--------','--------------','--------','--------------','--------','--------------','--------');

for i=1:numel(Nvec)
    N = Nvec(i);
    % For order 1
    L1h1 = L1h(1,i);
    if i>1, O1h = ord_h(1,i-1); else, O1h=NaN; end
    L1hu1 = L1hu(1,i);
    if i>1, O1hu=ord_hu(1,i-1); else, O1hu=NaN; end
    % For order 3
    L1h3 = L1h(2,i);
    if i>1, O3h = ord_h(2,i-1); else, O3h=NaN; end
    L1hu3 = L1hu(2,i);
    if i>1, O3hu = ord_hu(2,i-1); else, O3hu=NaN; end

    fprintf(fid,'%10d %14.4e %8.3f %14.4e %8.3f %14.4e %8.3f %14.4e %8.3f\n', ...
        N, L1h1, O1h, L1hu1, O1hu, L1h3, O3h, L1hu3, O3hu);
end

fclose(fid);

% Print to stdout
fprintf('Smooth accuracy test (N, L1(h), order, L1(hu), order):\n');
fprintf('N\t1_order_h\t1_ord\t1_order_hu\t1_ord\t3_order_h\t3_ord\t3_order_hu\t3_ord\n');
for i=1:numel(Nvec)
    N = Nvec(i);
    if i==1
        fprintf('%d\t%.4e\t NaN \t%.4e\t NaN \t%.4e\t NaN \t%.4e\t NaN\n',...
            N, L1h(1,i), L1hu(1,i), L1h(2,i), L1hu(2,i));
    else
        fprintf('%d\t%.4e\t%.3f\t%.4e\t%.3f\t%.4e\t%.3f\t%.4e\t%.3f\n',...
            N, L1h(1,i), ord_h(1,i-1), L1hu(1,i), ord_hu(1,i-1), L1h(2,i), ord_h(2,i-1), L1hu(2,i), ord_hu(2,i-1));
    end
end

% Plot log-log L1(h) vs dx with slope lines
fig = figure('Color','w','Position',[100 100 900 600]);
c = lines(7);
dxvec = 1./Nvec;

loglog(dxvec, L1h(1,:), 'o-', 'color', c(1,:), 'DisplayName', 'Order 1 Euler');
hold on
loglog(dxvec, L1h(2,:), 's-', 'color', c(2,:), 'DisplayName', 'Order 3 SSPRK3');

% Reference slopes lines with scale factor for visibility
refx = [dxvec(1) dxvec(end)];
refy1 = L1h(1,1)*(refx/refx(1)).^1;
refy3 = L1h(2,1)*(refx/refx(1)).^3;

loglog(refx, refy1, 'k--', 'LineWidth',1.5, 'DisplayName', 'Slope 1');
loglog(refx, refy3, 'k-.', 'LineWidth',1.5, 'DisplayName', 'Slope 3');

grid on
xlabel('dx (m)');
ylabel('L1 error in h');
title('Smooth accuracy test: L1 error in h vs dx with slope refs');
legend('Location','northwest','NumColumns',1);
set(gca,'FontSize',11)
hold off

exportgraphics(fig, fullfile(imgdir,'p3_smooth_accuracy.png'), 'Resolution',150)

%% 2. Numerical diffusion: pulse propagation, bed flat, periodic, tend=0.3, N=200
N = 200; dx = 1/N; bc = 'periodic'; tend = 0.3;
bedflat = make_bed('flat',[]);

x = linspace(dx/2,1-dx/2,N);  % cell-centres [0,1]

% Initial condition: eta=1+0.2*exp(-((x-0.5)/0.05)^2), u=0 => h=eta (flat bed), hu=0
initpulse = @(x,Bc) deal(1+0.2*exp(-((x-0.5)/0.05).^2), zeros(size(x)));

% Reference: order 3, ssprk3, cfl=0.2
cfgref2.N = Nref; cfgref2.xl = 0; cfgref2.xr = 1; cfgref2.bc = bc;
cfgref2.bed = bedflat; cfgref2.tend = tend;
cfgref2.scheme = 'hr'; cfgref2.order = 3; cfgref2.rk = 'ssprk3'; cfgref2.cfl = 0.2;
cfgref2.init = initpulse;
cfgref2.tsave = tend;
outref2 = swe1d_run(cfgref2);
% Coarsen ref solution by averaging cells matching coarse grid
cellsz = Nref / N;
eta_ref_coarse = mean(reshape(outref2.H{1}+outref2.B, cellsz, N),1);
x_coarse = linspace(dx/2, 1 - dx/2, N); 

% Order 1, Euler, cfl 0.5
cfg1 = cfgref2; 
cfg1.N = N; 
cfg1.order = 1;
cfg1.rk = 'euler';
cfg1.cfl = 0.5;
out1 = swe1d_run(cfg1);

% Order 3, SSPRK3, cfl 0.2
cfg3 = cfgref2;
cfg3.N = N;
cfg3.order = 3;
cfg3.rk = 'ssprk3';
cfg3.cfl = 0.2;
out3 = swe1d_run(cfg3);

% Extract eta for the pulse (eta=h+B)
eta1 = out1.H{1} + out1.B;
eta3 = out3.H{1} + out3.B;
% Use coarsened eta_ref with coarse x grid
eta_ref = eta_ref_coarse;

% Plot eta pulses
fig = figure('Color','w','Position',[100 100 900 600]);
plot(out1.x, eta1, '-.', 'LineWidth',1.5, 'Color', c(1,:), 'DisplayName', 'Order 1 Euler');
hold on
plot(out3.x, eta3, '-', 'LineWidth',1.5, 'Color', c(2,:), 'DisplayName', 'Order 3 SSPRK3');
plot(x_coarse, eta_ref, ':', 'LineWidth',2, 'Color', [0 0 0], 'DisplayName', 'Reference');

xlabel('x (m)')
ylabel('\eta (m)')
title('Numerical diffusion of pulse (bed flat, periodic)')
legend('Location','best')
grid on
set(gca,'FontSize',11)
hold off
exportgraphics(fig, fullfile(imgdir,'p3_pulse_diffusion.png'), 'Resolution',150)

% Peak amplitude (max eta-1) comparison
peak_ref = max(eta_ref - 1);
peak_1 = max(eta1 - 1);
peak_3 = max(eta3 - 1);

% Save peak to text file
fname2 = fullfile(resdir,'p3_pulse_peak.txt');
fid = fopen(fname2,'w');
fprintf(fid,'Pulse peak amplitude (max eta-1) at t=%.3f, N=%d\n',tend,N);
fprintf(fid,'Reference (order 3 fine grid): %.6f\n',peak_ref);
fprintf(fid,'Order 1 Euler: %.6f\n',peak_1);
fprintf(fid,'Order 3 SSPRK3: %.6f\n',peak_3);
fclose(fid);

fprintf('\nPulse amplitude peaks:\n  Reference = %.6f\n  Order 1 Euler = %.6f\n  Order 3 SSPRK3 = %.6f\n', peak_ref, peak_1, peak_3);

%% 3. Well-balancing: lake at rest over Gaussian bed, N=100, tend=2
Nwb = 100; tendwb = 2;
bedg = make_bed('gaussian',[0.4 0.5 0.05]);
bcwb = 'wall';

eta0 = 1; % lake at rest eta
init_lake = @(x,Bc) deal(max(0,eta0 - Bc), zeros(size(x)));

cfghr = struct;
cfghr.N = Nwb; cfghr.xl=0; cfghr.xr=1; cfghr.bc=bcwb; cfghr.bed = bedg;
cfghr.tend = tendwb; cfghr.scheme = 'hr'; cfghr.order = 3; cfghr.rk = 'ssprk3'; cfghr.cfl=0.2;
cfghr.init = init_lake; cfghr.tsave = tendwb;

out_hr = swe1d_run(cfghr);

cfgcentral = cfghr;
cfgcentral.scheme = 'central'; % naive, non well-balanced

out_central = swe1d_run(cfgcentral);

% max |u| of both solutions
u_hr = zeros(size(out_hr.H{1}));
wet_hr = out_hr.H{1} > 1e-10;
u_hr(wet_hr) = abs(out_hr.HU{1}(wet_hr)./out_hr.H{1}(wet_hr));
umax_hr = max(u_hr);

u_cent = zeros(size(out_central.H{1}));
wet_cent = out_central.H{1} > 1e-10;
u_cent(wet_cent) = abs(out_central.HU{1}(wet_cent)./out_central.H{1}(wet_cent));
umax_cent = max(u_cent);

% Save results
fname3 = fullfile(resdir,'p3_well_balancing.txt');
fid = fopen(fname3,'w');
fprintf(fid,'Well-balancing test, lake at rest, N=%d, tend=%.2f\n',Nwb,tendwb);
fprintf(fid,'Max |u| HR-3 scheme: %.6e\n', umax_hr);
fprintf(fid,'Max |u| Naive central-3 scheme: %.6e\n', umax_cent);
fclose(fid);

fprintf('\nWell-balancing max |u|:\n - HR-3: %.6e\n - Naive central-3: %.6e\n', umax_hr, umax_cent);

%% 4. Cost: print nsteps and nupdates for runs in 1 (N=200)
Ncost = 200;
cfgcost = struct;
cfgcost.N = Ncost;
cfgcost.xl = 0;
cfgcost.xr = 1;
cfgcost.bc = 'periodic';
cfgcost.bed = bed;
cfgcost.tend = tend;

% order1, euler, cfl=0.5
cfgcost1 = cfgcost;
cfgcost1.order = 1;
cfgcost1.rk = 'euler';
cfgcost1.cfl = 0.5;
cfgcost1.scheme = 'hr';
cfgcost1.init = initfun;
cfgcost1.tsave = tend;

outcost1 = swe1d_run(cfgcost1);

% order3, ssprk3, cfl=0.2
cfgcost3 = cfgcost;
cfgcost3.order = 3;
cfgcost3.rk = 'ssprk3';
cfgcost3.cfl = 0.2;
cfgcost3.scheme = 'hr';
cfgcost3.init = initfun;
cfgcost3.tsave = tend;

outcost3 = swe1d_run(cfgcost3);

fprintf('\nCost at N=200:\n');
fprintf('Order 1 Euler: nsteps = %d, nupdates = %d\n', outcost1.nsteps, outcost1.nupdates);
fprintf('Order 3 SSPRK3: nsteps = %d, nupdates = %d\n', outcost3.nsteps, outcost3.nupdates);

% Save cost to file
fname4 = fullfile(resdir,'p3_cost.txt');
fid = fopen(fname4,'w');
fprintf(fid,'Cost at N=200:\n');
fprintf(fid,'Order 1 Euler: nsteps = %d, nupdates = %d\n', outcost1.nsteps, outcost1.nupdates);
fprintf(fid,'Order 3 SSPRK3: nsteps = %d, nupdates = %d\n', outcost3.nsteps, outcost3.nupdates);
fclose(fid);

fprintf('\nEXP p3 DONE\n');
end
