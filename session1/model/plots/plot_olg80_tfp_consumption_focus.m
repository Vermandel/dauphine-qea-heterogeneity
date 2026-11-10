% Consumption on impact across ages, after the TFP shock (Session 1 slides).
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
% Standalone run: solve the OLG economy and the RANK benchmark first.
dynare olg80 -DSHOCK_TYPE=1 noclearall	% writes results_tfp.mat
dynare rank80 -DSHOCK_TYPE=1 noclearall	% writes results_rank_tfp.mat
% Then load the saved TFP experiments: OLG transition and RANK benchmark.
O = load('results_tfp.mat');
R = load('results_rank_tfp.mat');
oo = O.oo_;
on = cellstr(O.M_.endo_names);
rn = cellstr(R.M_.endo_names);
io = @(n) find(strcmp(on,n),1);
ir = @(n) find(strcmp(rn,n),1);
J = 80; age = (20:19+J)';
% Percentage deviation of each age consumption on impact, against the RANK response.
pdev = @(o,i,n) 100*(o.endo_simul(i,n)/o.steady_state(i)-1);
imp = arrayfun(@(j) pdev(oo,io(sprintf('c_%d',j)),2), 1:J);
rc  = pdev(R.oo_,ir('C'),2);
blue = [0 0.30 0.60]; red = [0.80 0.15 0.10];

figure('Color','w','Position',[100 100 520 640]);
bar(age,imp,1,'FaceColor',blue,'EdgeColor','none');
hold('on');
yline(rc,'--','Color',red,'LineWidth',2);
xline(64.5,':','retirement');
grid('on');
xlabel('Age');
ylabel('% deviation on impact');
title('Consumption on impact, by age');
legend({'OLG by age','RANK'},'Location','southwest');

% Export next to the model in the release layout, or to the session figures folder in the source tree.
if isfolder(fullfile(modeldir,'..','figures'))
	figdir = fullfile(modeldir,'..','figures');
else
	figdir = modeldir;
end
exportgraphics(gcf,fullfile(figdir,'olg80_tfp_consumption_focus.png'),'Resolution',300);
