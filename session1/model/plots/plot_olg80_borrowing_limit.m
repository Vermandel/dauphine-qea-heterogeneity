% Plot borrowing-limit effects on assets, MPCs, and transition dynamics.
% Run this file after adding Dynare 7.x to the MATLAB path.
% Locate the model files in the same folder (release download) or its parent.
here = fileparts(mfilename('fullpath'));
if isfile(fullfile(here,'olg80.mod'))
	modeldir = here;
	figdir = fullfile(modeldir,'figures'); % standalone model download
elseif isfile(fullfile(here,'..','olg80.mod'))
	modeldir = fullfile(here,'..');
	figdir = fullfile(modeldir,'..','figures'); % session1/figures, read by TeX
else
	error('Place this script next to the Session 1 model files or in model/plots.');
end
cd(modeldir);
dynare olg80 -DSHOCK_TYPE=1 noclearall      % results_tfp.mat
dynare rank80 -DSHOCK_TYPE=1 noclearall     % results_rank_tfp.mat
% olg80.mod writes the analytic unconstrained MPC to mpc_bc0.mat.
dynare olg80_bc noclearall                  % constrained mpc_bc1.mat
global M_ oo_
if ~isfolder(figdir)
	mkdir(figdir);
end
names = cellstr(M_.endo_names); pn = cellstr(M_.param_names);
J = 80; age = (20:19+J)'; t = 0:size(oo_.endo_simul,2)-1;
blue = [0 0.30 0.60]; orange = [0.85 0.45 0.05]; red = [0.80 0.15 0.10]; grey = [0.45 0.45 0.45 0.22];
z = arrayfun(@(j) M_.params(strcmp(pn,sprintf('z_%d',j))), 1:J)';
pdev = @(o,i) 100*(o.endo_simul(i,:)/o.steady_state(i)-1);

% Both MPCs hold prices fixed: free annuity formula versus constrained remaining-life budgets.
M0 = load('mpc_bc0.mat'); M1 = load('mpc_bc1.mat');
aprof = [M0.a_ss M1.a_ss]; mpc = [M0.mpc M1.mpc];
figure('Color','w','Position',[50 50 1100 420]); tiledlayout(1,2,'TileSpacing','compact');
nexttile; plot(age,aprof(:,1),'LineWidth',1.8); hold('on'); plot(age,aprof(:,2),'LineWidth',1.8);
yline(0,':k'); grid('on'); xlabel('Age'); ylabel('End-of-period assets'); title('Stationary asset profiles');
legend('No borrowing limit','Limit: a >= 0','Location','northwest');
nexttile; plot(age,mpc(:,1),'LineWidth',1.8); hold('on'); plot(age,mpc(:,2),'LineWidth',1.8);
grid('on'); ylim([0 1.05]); xlabel('Age'); ylabel('MPC');
title('Impact MPC out of a windfall');
legend('No borrowing limit','Limit: a >= 0','Location','north');
exportgraphics(gcf,fullfile(figdir,'olg80_borrowing_mpc.pdf'),'ContentType','vector');
fprintf('Binding ages: %s\n',mat2str(age(aprof(1:J,2)<1e-8 & age<age(end))'));
fprintf('MPC at ages 25, 35, 75, 99 (no limit): %s\n',mat2str(mpc([6 16 56 80],1)',6));
fprintf('MPC at ages 25, 35, 75, 99 (limit):    %s\n',mat2str(mpc([6 16 56 80],2)',6));

% TFP transition with and without the limit, against RANK.
O = load('results_tfp.mat'); R = load('results_rank_tfp.mat');
on = cellstr(O.M_.endo_names); rn = cellstr(R.M_.endo_names);
io = @(n) find(strcmp(on,n),1); ir = @(n) find(strcmp(rn,n),1); ib = @(n) find(strcmp(names,n),1);
assert(min(arrayfun(@(j) min(oo_.endo_simul(ib(sprintf('a_%d',j)),:)),1:J-1)) > -1e-8, 'Borrowing limit violated');
figure('Color','w','Position',[100 100 1000 600]); tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
econ = {O.oo_, oo_}; idx = {io, ib}; col = {blue, orange};
tO = 0:size(O.oo_.endo_simul,2)-1; tR = 0:size(R.oo_.endo_simul,2)-1; tq = {tO, t};
ttl = {'Consumption, no borrowing limit','Consumption, a_j \geq 0'};
lab = {'OLG aggregate','OLG aggregate, a_j \geq 0'};
for q = 1:2
	nexttile; hold('on');
	for j = 1:J, plot(tq{q},pdev(econ{q},idx{q}(sprintf('c_%d',j))),'Color',grey,'LineWidth',0.8); end
	h1 = plot(tq{q},pdev(econ{q},idx{q}('C')),'Color',col{q},'LineWidth',2.5);
	h2 = plot(tR,pdev(R.oo_,ir('C')),'--','Color',red,'LineWidth',2);
	h0 = plot(nan,nan,'Color',grey(1:3),'LineWidth',0.8);
	title(ttl{q}); ylabel('% deviation'); grid('on'); xlim([0 40]); ylim([-1.2 0.3]);
	legend([h0 h1 h2],{'OLG by age',lab{q},'RANK'},'Location','southeast');
end
vars = {'y','k','r','w'}; tt = {'Output','Capital','Net interest rate','Wage'};
for q = 1:4
	nexttile; hold('on');
	if q == 3, f = @(o,i) 100*(o.endo_simul(i,:)-o.steady_state(i)); yl = 'percentage points';
	else, f = pdev; yl = '% deviation'; end
	plot(tO,f(O.oo_,io(vars{q})),'Color',blue,'LineWidth',2.5);
	plot(t,f(oo_,ib(vars{q})),'Color',orange,'LineWidth',2.5);
	plot(tR,f(R.oo_,ir(vars{q})),'--','Color',red,'LineWidth',2);
	title(tt{q}); ylabel(yl); grid('on'); xlim([0 40]);
	if q == 1, legend('OLG','OLG, a_j \geq 0','RANK','Location','southeast'); end
end
xlabel(tl,'Years (unexpected shock at date 1)');
exportgraphics(gcf,fullfile(figdir,'olg80_tfp_transition_bc.pdf'),'ContentType','vector');
