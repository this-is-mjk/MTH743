function exp_pert
% Small perturbation lake at rest tests (LeVeque 1998) eps=0.2,0.001,
% N=100 methods, reference hr3 N=3200, all saved figures, tables, MP4

set(0,'DefaultFigureVisible','off');
s=settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue='light';

root = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

% Common params
cfg0.g = 9.81;
cfg0.xl = 0; cfg0.xr = 1;
cfg0.bc = 'open';
cfg0.tend = 0.2;

% Bed cosbump params from task
a = 0.25; x0 = 0.5; w = 0.1; 
bed = make_bed('cosbump',[a x0 w]);

% Grid resolutions
Ncoarse = 100;
Nref = 3200;
xcoarse = linspace(cfg0.xl+(0.5/Ncoarse),cfg0.xr-(0.5/Ncoarse),Ncoarse);
dx_ref = (cfg0.xr - cfg0.xl) / Nref;
xref = cfg0.xl + ((1:Nref) - 0.5)*dx_ref;

% Perturbation intervals
pert_interval = [0.1 0.2];

% Methods to test on coarse N, with schemes/orders/rk/cfl as specified
methods = {...
    struct('name','central-1','scheme','central','order',1,'rk','euler','cfl',0.5), ...
    struct('name','hr-1',     'scheme','hr',     'order',1,'rk','euler','cfl',0.5), ...
    struct('name','hr-3',     'scheme','hr',     'order',3,'rk','ssprk3','cfl',0.4)};

epsvals = [0.2,0.001];

% Reference solution config for each eps (same except N and CFL)
cfgref_base = struct('N',Nref,'xl',cfg0.xl,'xr',cfg0.xr,'bed',bed,'bc',cfg0.bc,...
    'scheme','hr','order',3,'rk','ssprk3','cfl',0.2,'tend',cfg0.tend);

% Output containers
out.ref = cell(2,1);
out.coarse = cell(length(epsvals),length(methods));
coords = xcoarse; % all coarse runs share xcoarse (1xNcoarse)

% Define line styles for later use in plotting
lstyles = {'-','--','-.'};

%% Define eta0 init functions with perturbation at given eps
% Lake at rest free surface eta = 1 except 1+eps over perturbation interval

init_eta_eps = @(eps) @(x,Bc) ...
    deal(max(0,(eta0_eps(x,eps) - Bc)), zeros(size(x)));

    function eta = eta0_eps(x,eps)
        eta = ones(size(x));
        idx = x > pert_interval(1) & x < pert_interval(2);
        eta(idx) = 1 + eps;
    end

%% Compute reference solutions and coarse solutions

fprintf('Starting reference solutions (N=3200)...\n');
for i=1:2
    e = epsvals(i);
    % Init with perturbation
    cfgref = cfgref_base;
    cfgref.init = init_eta_eps(e);
    % Use high-resolution grid
    cfgref.N = Nref;
    out.ref{i} = swe1d_run(cfgref);
    % Reference height at time tend (coarse averages will be done later)
end

fprintf('Reference runs done.\n');

% Cell average reference on coarse grid by reshaping + mean in blocks
avg_ref_eta_scaled = cell(2,1);
for i=1:2
    Hr = out.ref{i}.H{end};
    % average to coarse grid by reshaping blocks of (Nref/Ncoarse)=32
    blocksize = Nref/Ncoarse;
    % cell-averaged SURFACE: block mean of (h + B) on the fine grid. Averaging only h and adding the
    % bed at the coarse cell centre would leave an error mean(B)-B(xc) ~ B''dx^2/24 (~5e-4 on the cosine bump),
    % i.e. half of the signal for eps = 0.001.
    Bfine = bed.f(out.ref{i}.x,0);
    eta_ref = mean(reshape(Hr + Bfine,blocksize,Ncoarse),1);
    % scale (eta-1)/eps
    ep = epsvals(i);
    avg_ref_eta_scaled{i} = (eta_ref - 1)/ep;
end

%% Run coarse runs for each eps, for each method
for i=1:2
    ep = epsvals(i);
    fprintf('Running coarse N=100 for eps=%.4f\n',ep);
    for m=1:length(methods)
        meth = methods{m};
        cfg = struct;
        cfg.N = Ncoarse;
        cfg.xl = cfg0.xl; cfg.xr = cfg0.xr;
        cfg.bed = bed;
        cfg.bc = cfg0.bc;
        cfg.g = cfg0.g;
        cfg.scheme = meth.scheme;
        cfg.order = meth.order;
        cfg.rk = meth.rk;
        cfg.cfl = meth.cfl;
        cfg.tend = cfg0.tend;
        cfg.init = init_eta_eps(ep);
        outm = swe1d_run(cfg);
        out.coarse{i,m} = outm;
    end
end

%% 1. Figures for each eps: plot (eta-1)/eps at tend, 3 methods + ref + bed scaled

colsline = lines(7);
grey = 0.5*ones(1,3);
figh = cell(2,1);

for i=1:2
    ep = epsvals(i);
    fig = figure('Color','w','Position',[100 100 900 600]);
    hold on; box on;
    % plot methods scaled eta
    for m=1:length(methods)
        outm = out.coarse{i,m};
        Hf = outm.H{end};
        Bcoarse = bed.f(coords,0);
        eta = Hf + Bcoarse;
        eta_scaled = (eta - 1)/ep;
        plot(coords,eta_scaled,'LineWidth',1.8,'Color',colsline(m,:),...
            'LineStyle',lstyles{m},'DisplayName',methods{m}.name);
    end
    % plot ref scaled (black line)
    plot(coords, avg_ref_eta_scaled{i}, 'k-', 'LineWidth',2, 'DisplayName','reference');
    % plot scaled bed (grey)
    Bcoarse = bed.f(coords,0);
    bed_scaled = (Bcoarse - min(Bcoarse))/max(Bcoarse - min(Bcoarse)); % normalize
    bed_scaled = 0.25*bed_scaled; % scale bed for figure
    plot(coords, bed_scaled, '-', 'Color', grey, 'LineWidth', 1.5, 'DisplayName','bed (scaled)');
    xlabel('x [m]','FontSize',13);
    ylabel('(η - 1) / ε','FontSize',13);
    title(sprintf('Scaled free surface perturbation at t=%.2f, \\varepsilon=%.4g',cfg0.tend,ep),'FontSize',14);
    legend('Location','northwest','FontSize',12,'NumColumns',1);
    set(gca,'FontSize',12);
    ylim([-0.05 2.1]);
    figh{i} = fig;
    fname = fullfile(imgdir,sprintf('pert_eps%04d.png',round(ep*1e4)));
    exportgraphics(fig,fname,'Resolution',150);
    close(fig);
end

%% 2. Table: eps, method, L1 error and max|(eta-1)/eps|, and ratio for eps=0.001 (central-1 error)/(hr-3 error)

fid = fopen(fullfile(resdir,'pert.txt'),'w');
fprintf(fid,'eps\tmethod\tL1_error\tmax_abs_scaled_eta\n');

errmat = nan(2,length(methods));
maxmat = nan(2,length(methods));

% Get method names safely (fix of dot indexing on cell)
mnames = cellfun(@(m) m.name, methods, 'UniformOutput', false);

for i=1:2
    ep = epsvals(i);
    refv = avg_ref_eta_scaled{i};
    dx = (cfg0.xr - cfg0.xl)/Ncoarse;
    for m=1:length(methods)
        outm = out.coarse{i,m};
        Hf = outm.H{end};
        Bcoarse = bed.f(coords,0);
        eta = Hf + Bcoarse;
        eta_scaled = (eta - 1)/ep;
        err = sum(abs(eta_scaled - refv)) * dx;
        mx = max(abs(eta_scaled));
        errmat(i,m) = err;
        maxmat(i,m) = mx;
        fprintf(fid,'%.4g\t%s\t%.4e\t%.5f\n', ep, mnames{m}, err, mx);
    end
end

% Ratio of errors for eps=0.001 case (index 2)
idx_eps1 = 2;
idx_central = find(strcmp(mnames,'central-1'));
idx_hr3 = find(strcmp(mnames,'hr-3'));
if ~isempty(idx_central) && ~isempty(idx_hr3)
    ratio_err = (errmat(idx_eps1,idx_central)-errmat(idx_eps1,idx_hr3))/errmat(idx_eps1,idx_hr3);
    fprintf(fid,'\nError ratio (central-1 - hr-3)/hr-3 for eps=0.001: %.4f\n', ratio_err);
    fprintf('Error ratio (central-1 - hr-3)/hr-3 for eps=0.001: %.4f\n', ratio_err);
else
    ratio_err = NaN;
end

fclose(fid);

%% 3. Convergence of hr-3 for N=50,100,200,400 for eps=0.001 

Nlist = [50 100 200 400];
ep = 0.001;
order_est = nan(1,length(Nlist)-1);
L1errs = nan(1,length(Nlist));

fprintf('Convergence study for hr-3, eps=0.001\n');
for k=1:length(Nlist)
    Nk = Nlist(k);
    cfg = struct;
    cfg.N = Nk;
    cfg.xl = cfg0.xl; cfg.xr = cfg0.xr;
    cfg.bed = bed;
    cfg.bc = cfg0.bc;
    cfg.g = cfg0.g;
    cfg.scheme = 'hr';
    cfg.order = 3;
    cfg.rk = 'ssprk3';
    cfg.cfl = 0.4;
    cfg.tend = cfg0.tend;
    cfg.init = init_eta_eps(ep);
    outk = swe1d_run(cfg);
    % Reference data downsampled or interpolated for error
    % For consistency, use averaging from ref N=3200 to Nk
    nxref = Nref;
    blocksize = nxref / Nk;
    Hr_ref = out.ref{2}.H{end};
    xk = cfg.xl + ((1:Nk)-0.5)*( (cfg.xr-cfg.xl)/Nk );
    Bk = bed.f(xk,0);
    Bfine = bed.f(out.ref{2}.x,0);
    eta_ref_k = mean(reshape(Hr_ref + Bfine,blocksize,Nk),1);   % block mean of the surface
    Hnum = outk.H{end};
    eta_num = Hnum + Bk;
    eta_scaled_ref = (eta_ref_k - 1)/ep;
    eta_scaled_num = (eta_num - 1)/ep;
    dxk = (cfg.xr-cfg.xl)/Nk;
    errk = sum(abs(eta_scaled_num - eta_scaled_ref))*dxk;
    L1errs(k) = errk;
    fprintf('N=%4d  L1 error=%.4e\n', Nk, errk);
end

% Compute order of convergence log2(eN/e2N)
for k=1:length(Nlist)-1
    order_est(k) = log2(L1errs(k)/L1errs(k+1));
end

fid = fopen(fullfile(resdir,'pert_convergence.txt'),'w');
fprintf(fid,'N\tL1_error\tObserved_order\n');
fprintf(fid,'%d\t%.6e\t-\n', Nlist(1), L1errs(1));
for k=2:length(Nlist)
    fprintf(fid,'%d\t%.6e\t%.4f\n', Nlist(k), L1errs(k), order_est(k-1));
end
fclose(fid);

%% 4. MP4 pert_eps0001.mp4 eps=0.001, N=200, hr-3 vs central-1, two panels, 15fps, 81 frames tsave linspace(0,0.2,81)

ep = 0.001; Nmp = 200; nts = 81;
tlist = linspace(0,cfg0.tend,nts);

cfghr = struct('N',Nmp,'xl',cfg0.xl,'xr',cfg0.xr,'bed',bed,'bc',cfg0.bc,'g',cfg0.g,...
    'scheme','hr','order',3,'rk','ssprk3','cfl',0.4,'tend',cfg0.tend,'tsave',tlist,'init',init_eta_eps(ep));
out_hr = swe1d_run(cfghr);

cfgcen = struct('N',Nmp,'xl',cfg0.xl,'xr',cfg0.xr,'bed',bed,'bc',cfg0.bc,'g',cfg0.g,...
    'scheme','central','order',1,'rk','euler','cfl',0.5,'tend',cfg0.tend,'tsave',tlist,'init',init_eta_eps(ep));
out_cen = swe1d_run(cfgcen);

xmp = out_hr.x;
Bmp = bed.f(xmp,0);

% Set up video
vidname = fullfile(root,'outputs','figures','pert_eps0001.mp4');
v = VideoWriter(vidname,'MPEG-4');
v.FrameRate = 15;
open(v);

colsvid = lines(3);
grey = 0.5*ones(1,3);

figv = figure('Color','w','Position',[100 100 900 600]);

for k=1:nts
    clf;
    % Two panels
    % left: hr-3
    subplot(1,2,1);
    eta_hr = out_hr.H{k} + Bmp;
    plot(xmp, (eta_hr - 1)/ep, 'Color', colsvid(1,:), 'LineWidth', 1.8);
    hold on;
    bed_scaledv = (Bmp - min(Bmp))/max(Bmp - min(Bmp)); % scale bed
    bed_scaledv = 0.3*bed_scaledv;
    plot(xmp, bed_scaledv, 'Color', grey, 'LineWidth', 1.4);
    ylim([-0.05 2.1]);
    box on;
    xlabel('x [m]');
    ylabel('(η - 1)/ε');
    title(sprintf('hr-3 solution at t=%.3f', tlist(k)));
    set(gca,'FontSize',12);

    % right: central-1
    subplot(1,2,2);
    eta_cen = out_cen.H{k} + Bmp;
    plot(xmp,(eta_cen - 1)/ep,'Color',colsvid(2,:),'LineWidth',1.8);
    hold on;
    plot(xmp, bed_scaledv, 'Color', grey, 'LineWidth', 1.4);
    ylim([-0.05 2.1]);
    box on;
    xlabel('x [m]');
    ylabel('(η - 1)/ε');
    title(sprintf('central-1 solution at t=%.3f', tlist(k)));
    set(gca,'FontSize',12);

    drawnow;
    frame = getframe(figv);
    writeVideo(v,frame);
end
close(v);
close(figv);

fprintf('EXP pert DONE\n');

end
