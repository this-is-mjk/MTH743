% exp_beds: reference images (no simulation)
%   beds_reference.png  - every bed profile used by the experiments, with the initial water surface,
%                         the experiment that uses it and the steepest slope max|B'|
%   hr_schematic.png    - what the hydrostatic reconstruction does at one cell interface
set(0,'DefaultFigureVisible','off');
s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light';
root   = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

colBed = [0.55 0.38 0.2]; colWater = [0.62 0.80 0.95]; colSurf = [0 0.447 0.741];

% ---- 1. gallery of the beds ---------------------------------------------------------------------------
% name, bed, domain, initial surface eta0(x) (function handle of x and B), where used
E = {
 'p1, p2: Gaussian [a x0 s] = [0.4 5 0.5]', make_bed('gaussian',[0.4 5 0.5]), [0 10], @(x,B) ones(size(x)),               'lake at rest, eta = 1';
 'p2: cosine bump [a x0 w] = [0.5 5 2]',                      make_bed('cosbump', [0.5 5 2]),   [0 10], @(x,B) ones(size(x)),               'lake at rest, eta = 1';
 'p2: tall narrow Gaussian [0.8 3 0.4]',                      make_bed('gaussian',[0.8 3 0.4]), [0 10], @(x,B) ones(size(x)),               'lake at rest, eta = 1';
 'pert, rk: cosine bump [0.5 0.5 0.1]',                       make_bed('cosbump', [0.5 0.5 0.1]), [0 1], @(x,B) 1+0.2*(x>0.1 & x<0.2),      'eta = 1 + eps on 0.1<x<0.2 (eps = 0.2 shown)';
 'p3 accuracy: sin^2(pi x)',                                   make_bed('sin2',[]),               [0 1], @(x,B) B + 5 + exp(cos(2*pi*x)),    'smooth flow, h = 5 + exp(cos 2 pi x)';
 'p3 well-balancing: narrow Gaussian [0.4 0.5 0.05]',         make_bed('gaussian',[0.4 0.5 0.05]), [0 1], @(x,B) ones(size(x)),             'lake at rest, eta = 1';
 'dam: step [B0 B1 xs] = [0 0.5 12]',                          make_bed('step',[0 0.5 12]),      [0 20], @(x,B) B + 2*(x<8) + 1*(x>=8 & x<12) + 0.5*(x>=12), 'depths 2 | 1 | 0.5 (surface 2 | 1 | 1)';
 'tsunami: shelf [x1 x2 b1 b2] = [30 90 -10 -1]',             make_bed('slope',[30 90 -10 -1]), [0 100], @(x,B) 0.05*exp(-((x-15)/3).^2),   'eta = 0.05 exp(-((x-15)/3)^2) (5 cm hump)';
};
fig = figure('Color','w','Position',[50 50 1500 1000]);
tl = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
fid = fopen(fullfile(resdir,'beds.txt'),'w');
fprintf(fid,'bed\tdomain length L\ttypical depth H\tL/H\tmax |B''|\n');
for k = 1:size(E,1)
  nexttile; hold on; box on;
  bed = E{k,2}; dom = E{k,3}; x = linspace(dom(1),dom(2),2000); B = bed.f(x,0); eta = E{k,4}(x,B);
  ymin = min(B) - 0.15*(max(eta)-min(B));
  fill([x fliplr(x)],[B ymin*ones(size(B))],colBed,'EdgeColor','none');
  fill([x fliplr(x)],[eta fliplr(B)],colWater,'EdgeColor','none','FaceAlpha',0.6);
  plot(x,B,'Color',[0.3 0.2 0.1],'LineWidth',1.8);
  plot(x,eta,'Color',colSurf,'LineWidth',1.6);
  slope = max(abs(diff(B)./diff(x)));
  H = mean(eta - B); L = dom(2)-dom(1);
  xlim(dom); ylim([ymin, max(eta)+0.12*(max(eta)-ymin)]);
  xlabel('x [m]'); ylabel('elevation [m]'); set(gca,'FontSize',11);
  title({E{k,1}, sprintf('%s;  max|B''| = %.2f', E{k,5}, slope)},'FontSize',10,'Interpreter','none');
  fprintf(fid,'%s\t%g\t%.3g\t%.3g\t%.3g\n',E{k,1},L,H,L/H,slope);
end
% 2D bump
nexttile; xg = linspace(0,1,200); [X,Y] = meshgrid(xg,xg); B2 = 0.4*exp(-((X-0.5).^2+(Y-0.5).^2)/0.02);
contourf(X,Y,B2,12,'LineColor','none'); colormap(gca,flipud(gray)); cb = colorbar; cb.Label.String = 'B [m]'; axis equal tight;
xlabel('x [m]'); ylabel('y [m]'); set(gca,'FontSize',11);
[gx,gy] = gradient(B2,xg(2)-xg(1)); sl2 = max(hypot(gx(:),gy(:)));
title({'2D: B = 0.4 exp(-r^2/0.02), r from (0.5,0.5)', sprintf('lake at rest, eta = 1;  max|grad B| = %.2f',sl2)},'FontSize',10);
fprintf(fid,'2D Gaussian bump\t1\t~0.7\t~1.4\t%.3g\n',sl2);
fclose(fid);
title(tl,'Beds used in the experiments (brown = bed, blue = water, line = free surface eta)','FontSize',13);
exportgraphics(fig, fullfile(imgdir,'beds_reference.png'),'Resolution',150); close(fig);

% ---- 2. hydrostatic reconstruction schematic -------------------------------------------------------
fig = figure('Color','w','Position',[50 50 1400 560]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
cases = { struct('title','(a) lake at rest: eta_L = eta_R = 1','etaL',1.0,'etaR',1.0), ...
          struct('title','(b) a real step in the surface: eta_L = 1.2, eta_R = 0.9','etaL',1.2,'etaR',0.9) };
BL = 0; BR = 0.4; Bs = max(BL,BR);
for c = 1:2
  cs = cases{c}; nexttile; hold on; box on;
  fill([0 1 1 0],[-0.15 -0.15 BL BL],colBed,'EdgeColor','none');          % left bed
  fill([1 2 2 1],[-0.15 -0.15 BR BR],colBed,'EdgeColor','none');          % right bed
  fill([0 1 1 0],[BL BL cs.etaL cs.etaL],colWater,'EdgeColor','none');    % left water (true state)
  fill([1 2 2 1],[BR BR cs.etaR cs.etaR],colWater,'EdgeColor','none');    % right water (true state)
  hLs = max(0,cs.etaL-Bs); hRs = max(0,cs.etaR-Bs);
  % reconstructed states: water above the common bed B*
  fill([0.42 0.98 0.98 0.42],[Bs Bs cs.etaL cs.etaL],[0.20 0.55 0.85],'EdgeColor','none','FaceAlpha',0.55);
  fill([1.02 1.58 1.58 1.02],[Bs Bs cs.etaR cs.etaR],[0.20 0.55 0.85],'EdgeColor','none','FaceAlpha',0.55);
  plot([0 2],[Bs Bs],'k--','LineWidth',1.3);
  plot([1 1],[-0.15 1.45],'k:','LineWidth',1);
  text(1.08,Bs-0.07,'B^* = max(B_L, B_R) = B_R','FontSize',11,'Color','w','FontWeight','bold');
  text(0.05,-0.08,'bed B_L','FontSize',11,'Color','w'); text(1.05,-0.08,'bed B_R (higher)','FontSize',11,'Color','w');
  % depth brackets
  plot([0.30 0.30],[BL cs.etaL],'k-','LineWidth',1.5); text(0.05,(BL+cs.etaL)/2,sprintf('h_L = %.1f',cs.etaL-BL),'FontSize',11,'BackgroundColor','w');
  plot([0.36 0.36],[Bs cs.etaL],'-','Color',[0 0.3 0.6],'LineWidth',2.5); text(0.44,(Bs+cs.etaL)/2,sprintf('h_L^* = %.1f',hLs),'FontSize',11,'Color',[0 0.2 0.5],'FontWeight','bold');
  plot([1.70 1.70],[BR cs.etaR],'k-','LineWidth',1.5); text(1.73,(BR+cs.etaR)/2+0.12,sprintf('h_R = %.1f',cs.etaR-BR),'FontSize',11,'BackgroundColor','w');
  plot([1.64 1.64],[Bs cs.etaR],'-','Color',[0 0.3 0.6],'LineWidth',2.5); text(1.04,(Bs+cs.etaR)/2-0.1,sprintf('h_R^* = %.1f',hRs),'FontSize',11,'Color',[0 0.2 0.5],'FontWeight','bold');
  text(1.0,1.40,'interface','HorizontalAlignment','center','FontSize',11);
  xlim([0 2]); ylim([-0.15 2.15]); set(gca,'XTick',[],'FontSize',11); ylabel('elevation [m]');
  if c == 1
    str = {'reconstructed states are EQUAL (h_L^* = h_R^* = 0.6, u = 0)','=> HLL flux = (0, g h*^2/2): no numerical diffusion, no mass flux', ...
           'pressure corrections g/2 (h^2 - h*^2) restore the true pressure of each cell','=> the cell sees g h^2/2 on both faces: flux difference = 0 = source'};
  else
    str = {'reconstructed states differ (h_L^* = 0.8, h_R^* = 0.5): a genuine Riemann problem','=> HLL creates a wave (the surface step is a real dam break)', ...
           'the bed step only acts through the corrections g/2 (h^2 - h*^2)','(hydrostatic force of the bed step on the water)'};
  end
  text(0.03,2.10,str,'FontSize',10,'VerticalAlignment','top','BackgroundColor',[1 1 1],'EdgeColor',[0.7 0.7 0.7]);
  title(cs.title,'FontSize',12);
end
title(tl,'Hydrostatic reconstruction at one interface: light blue = true water, dark blue = water above the common bed B^*','FontSize',13);
exportgraphics(fig, fullfile(imgdir,'hr_schematic.png'),'Resolution',150); close(fig);
fprintf('EXP beds DONE\n');
