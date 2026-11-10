% Plot the OLG calibration, life-cycle profiles, and asset-market equilibrium.
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
dynare olg80 -DSHOCK_TYPE=1 noclearall
global M_ oo_
% Export next to the model in the release layout, or to the session figures folder in the source tree.
if isfolder(fullfile(modeldir,'..','figures'))
	figdir = fullfile(modeldir,'..','figures');
else
	figdir = modeldir;
end
names = cellstr(M_.endo_names); pn = cellstr(M_.param_names);
J = 80; age = (20:19+J)';
z  = arrayfun(@(j) M_.params(strcmp(pn,sprintf('z_%d',j))), 1:J)';
Nj = arrayfun(@(j) M_.params(strcmp(pn,sprintf('N_%d',j))), 1:J)';
ss = @(n) oo_.steady_state(strcmp(names,n));
css = arrayfun(@(j) ss(sprintf('c_%d',j)), 1:J)';
ass = [arrayfun(@(j) ss(sprintf('a_%d',j)), 1:J-1)'; 0];
blue = [0 0.30 0.60]; red = [0.80 0.15 0.10];

figure('Color','w','Position',[100 100 760 520]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(age,z,'Color',blue,'LineWidth',2); grid('on');
xlabel('Age'); title('Labor efficiency e_j'); xline(65,':','retirement');
nexttile; plot(age,Nj,'Color',blue,'LineWidth',2); grid('on');
xlabel('Age'); title('Cohort mass N_j'); ylim([0 2/J]);
nexttile; plot(age(1:J-1),ones(J-1,1),'Color',blue,'LineWidth',2); hold('on');
plot(age(J),0,'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',5); grid('on');
xlabel('Age'); title('Survival rate s_j'); ylim([0 1.1]);
nexttile; plot(age,ss('w')*z,'Color',blue,'LineWidth',2); hold('on');
plot(age,css,'--','Color',red,'LineWidth',2); grid('on'); xlabel('Age');
title('Labor income and consumption');
legend('Labor income w e_j','Consumption c_j','Location','east');
exportgraphics(gcf,fullfile(figdir,'olg80_calibration_profiles.pdf'),'ContentType','vector');

figure('Color','w','Position',[100 100 920 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(age,css,'LineWidth',2); grid('on'); xlabel('Age'); ylabel('Consumption');
title('Steady-state consumption profile');
nexttile; plot(age,ass,'LineWidth',2); hold('on'); yline(0,'k:'); grid('on');
xlabel('Age'); ylabel('Assets'); title('Steady-state wealth profile');
exportgraphics(gcf,fullfile(figdir,'olg80_lifecycle.pdf'),'ContentType','vector');

% Asset supply and capital demand as functions of the interest rate
Lb = M_.params(strcmp(pn,'Lbar'));		% aggregate efficient labor
al = M_.params(strcmp(pn,'alpha'));		% capital share
be = M_.params(strcmp(pn,'beta'));		% discount factor
de = M_.params(strcmp(pn,'delta'));		% depreciation rate
si = M_.params(strcmp(pn,'sigma'));		% inverse intertemporal elasticity

rr = linspace(-0.015,0.04,120);			% grid of net interest rates
A  = nan(size(rr));				% asset supply on the grid
Kd = A;						% capital demand on the grid
h  = (0:J-1)';					% years since entering the economy

for i = 1:numel(rr)
	% Firm side: capital demand and wage implied by the interest rate
	R     = 1+rr(i);
	Kd(i) = (al*Lb^(1-al)/(rr(i)+de))^(1/(1-al));
	w     = (1-al)*Kd(i)^al*Lb^(-al);

	% Household side: Euler growth factor and consumption from the lifetime budget
	g = (be*R)^(1/si);
	c = sum(w*z./R.^h)/sum(g.^h./R.^h)*g.^h;

	% Assets accumulate with the budget constraint
	a = zeros(J,1);
	a(1) = w*z(1)-c(1);
	for j = 2:J
		a(j) = R*a(j-1)+w*z(j)-c(j);
	end

	% Aggregate asset supply
	A(i) = Nj'*a;
end

figure('Color','w','Position',[100 100 620 430]); hold('on');
plot(A,100*rr,'Color',blue,'LineWidth',2.2); plot(Kd,100*rr,'k','LineWidth',1.6);
plot(ss('k'),100*ss('r'),'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',7);
grid('on'); xlabel('Assets / capital per capita'); ylabel('Net interest rate r (%)'); xlim([0 12]);
legend('Asset supply','Capital demand','Equilibrium','Location','southwest');
exportgraphics(gcf,fullfile(figdir,'olg80_asset_market.pdf'),'ContentType','vector');
