% Compare OLG and RANK transition dynamics after the Session 1 shocks.
% Run this file after adding Dynare 7.x to the MATLAB path.
% Locate the model files in the same folder (release download) or its parent.
here = fileparts(mfilename('fullpath'));
if isfile(fullfile(here,'olg80.mod'))
	modeldir = here;
elseif isfile(fullfile(here,'..','olg80.mod'))
	modeldir = fullfile(here,'..');
else
	error('Place this script next to the Session 1 model files or in model/plots.');
end
cd(modeldir);
%% SHOCK_TYPE = 1
cd(modeldir);
dynare olg80 -DSHOCK_TYPE=1 noclearall
dynare rank80 -DSHOCK_TYPE=1 noclearall
global M_ oo_
% Export next to the model in the release layout, or to the session figures folder in the source tree.
if isfolder(fullfile(modeldir,'..','figures'))
	figdir = fullfile(modeldir,'..','figures');
else
	figdir = modeldir;
end
    O = load('results_tfp.mat'); figname = 'olg80_tfp_transition.pdf';
oo = O.oo_; J = 80;
on = cellstr(O.M_.endo_names); rn = cellstr(M_.endo_names);
io = @(n) find(strcmp(on,n),1); ir = @(n) find(strcmp(rn,n),1);
pdev = @(o,i) 100*(o.endo_simul(i,:)/o.steady_state(i)-1);
ppt  = @(o,i) 100*(o.endo_simul(i,:)-o.steady_state(i));
t = 0:size(oo.endo_simul,2)-1;
grey = [0.45 0.45 0.45 0.22]; blue = [0 0.30 0.60]; red = [0.80 0.15 0.10];

figure('Color','w','Position',[100 100 1000 600]);
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile([2 1]);
hold('on');
for j = 1:J, plot(t,pdev(oo,io(sprintf('c_%d',j))),'Color',grey,'LineWidth',0.8); end
h1 = plot(t,pdev(oo,io('C')),'Color',blue,'LineWidth',2.5);
h2 = plot(t,pdev(oo_,ir('C')),'--','Color',red,'LineWidth',2);
h0 = plot(nan,nan,'Color',grey(1:3),'LineWidth',0.8);
title('Consumption'); ylabel('% deviation'); grid('on'); xlim([0 40]);
legend([h0 h1 h2],{'OLG by age','OLG aggregate','RANK'},'Location','southeast');
ylim([-0.8 0.3]);
vars = {'y','k','r','w'}; ttl = {'Output','Capital','Net interest rate','Wage'};
for q = 1:4
    nexttile; hold('on');
    if q == 3, f = ppt; yl = 'percentage points'; else, f = pdev; yl = '% deviation'; end
    plot(t,f(oo,io(vars{q})),'Color',blue,'LineWidth',2.5);
    plot(t,f(oo_,ir(vars{q})),'--','Color',red,'LineWidth',2);
    title(ttl{q}); ylabel(yl); grid('on'); xlim([0 40]);
    if q == 1, legend('OLG','RANK','Location','southeast'); end
end
xlabel(tl,'Years (unexpected shock at date 1)');
exportgraphics(gcf,fullfile(figdir,figname),'ContentType','vector');

%% SHOCK_TYPE = 2
cd(modeldir);
dynare olg80 -DSHOCK_TYPE=2 noclearall
dynare rank80 -DSHOCK_TYPE=2 noclearall
% Export next to the model in the release layout, or to the session figures folder in the source tree.
if isfolder(fullfile(modeldir,'..','figures'))
	figdir = fullfile(modeldir,'..','figures');
else
	figdir = modeldir;
end
    O = load('results_dep.mat'); figname = 'olg80_depreciation.pdf';
oo = O.oo_; J = 80; age = (20:19+J)';
on = cellstr(O.M_.endo_names); rn = cellstr(M_.endo_names);
io = @(n) find(strcmp(on,n),1); ir = @(n) find(strcmp(rn,n),1);
pdev = @(o,i) 100*(o.endo_simul(i,:)/o.steady_state(i)-1);
ppt  = @(o,i) 100*(o.endo_simul(i,:)-o.steady_state(i));
t = 0:size(oo.endo_simul,2)-1;
grey = [0.45 0.45 0.45 0.22]; blue = [0 0.30 0.60]; red = [0.80 0.15 0.10];

figure('Color','w','Position',[100 100 1000 600]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
hold('on');
for j = 1:J, plot(t,pdev(oo,io(sprintf('c_%d',j))),'Color',grey,'LineWidth',0.8); end
h1 = plot(t,pdev(oo,io('C')),'Color',blue,'LineWidth',2.5);
h2 = plot(t,pdev(oo_,ir('C')),'--','Color',red,'LineWidth',2);
h0 = plot(nan,nan,'Color',grey(1:3),'LineWidth',0.8);
title('Consumption'); ylabel('% deviation'); grid('on'); xlim([0 40]);
legend([h0 h1 h2],{'OLG by age','OLG aggregate','RANK'},'Location','southeast');
nexttile; hold('on');
plot(t,pdev(oo,io('k')),'Color',blue,'LineWidth',2.5); plot(t,pdev(oo_,ir('k')),'--','Color',red,'LineWidth',2);
title('Capital'); ylabel('% deviation'); grid('on'); xlim([0 40]); legend('OLG','RANK','Location','southeast');
nexttile; hold('on');
plot(t,ppt(oo,io('r')),'Color',blue,'LineWidth',2.5); plot(t,ppt(oo_,ir('r')),'--','Color',red,'LineWidth',2);
title('Net interest rate'); ylabel('percentage points'); grid('on'); xlim([0 40]);
nexttile; hold('on');
imp = arrayfun(@(j) 100*(oo.endo_simul(io(sprintf('c_%d',j)),2)/oo.steady_state(io(sprintf('c_%d',j)))-1),1:J);
bar(age,imp,1,'FaceColor',blue,'EdgeColor','none'); grid('on');
rc = pdev(oo_,ir('C')); yline(rc(2),'--','Color',red,'LineWidth',2);
xline(64.5,':','retirement'); xlabel('Age'); ylabel('% deviation on impact'); title('Consumption on impact, by age');
xlabel(tl,'Years (unexpected shock at date 1)');
exportgraphics(gcf,fullfile(figdir,figname),'ContentType','vector');
