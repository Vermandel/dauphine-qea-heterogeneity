% Dauphine PSL, Master QEA — Gauthier Vermandel.
% Session 1: run OLG and RANK, then compare their aggregate shock responses.
% Place this script beside olg80.mod, rank80.mod and olg80_interest_rate.m.
% Add Dynare 7.2 to the MATLAB path before running.
shock_type   = 1;	% 1: productivity shock; 2: depreciation shock.
rerun_models = true;	% Set false to compare two MAT files already saved here.
modeldir     = fileparts(mfilename('fullpath'));	% Folder containing the course models.
cd(modeldir);

% Select the two files written by the same shock experiment.
if shock_type == 1
	olg_file  = 'results_tfp.mat';	% OLG productivity response.
	rank_file = 'results_rank_tfp.mat';	% RANK productivity response.
elseif shock_type == 2
	olg_file  = 'results_dep.mat';	% OLG depreciation response.
	rank_file = 'results_rank_dep.mat';	% RANK depreciation response.
else
	error('Choose shock_type = 1 or 2.');
end

if rerun_models
	% OLG saves the reference path; RANK reads it and draws the comparison.
	if shock_type == 1
		dynare olg80 -DSHOCK_TYPE=1 noclearall
		dynare rank80 -DSHOCK_TYPE=1 noclearall
	else
		dynare olg80 -DSHOCK_TYPE=2 noclearall
		dynare rank80 -DSHOCK_TYPE=2 noclearall
	end
else
	% Reuse saved paths without solving either economy again.
	assert(isfile(olg_file) && isfile(rank_file), ...
		'Run once with rerun_models = true to create both MAT files.');
	O = load(olg_file);	% OLG transition and stationary equilibrium.
	A = load(rank_file);	% RANK transition and stationary equilibrium.
	oo = O.oo_;	% OLG results.
	ra = A.oo_;	% RANK results.
	on = cellstr(O.M_.endo_names);	% OLG endogenous variables.
	rn = cellstr(A.M_.endo_names);	% RANK endogenous variables.
	io = @(n) find(strcmp(on,n),1);	% Row of an OLG aggregate.
	ir = @(n) find(strcmp(rn,n),1);	% Row of a RANK aggregate.
	time = 0:size(oo.endo_simul,2)-1;	% Years since the shock.
	vars = {'C','y','k','r'};	% Consumption, output, capital, return.
	titles = {'Consumption','Output','Capital','Net interest rate'};	% Panel names.

	figure('Name','OLG vs RANK');
	for q = 1:numel(vars)
		subplot(2,2,q);
			i = io(vars{q});	% OLG row for this aggregate.
			j = ir(vars{q});	% RANK row for this aggregate.
			if strcmp(vars{q},'r')
				x = 100*(oo.endo_simul(i,:)-oo.steady_state(i));
				y = 100*(ra.endo_simul(j,:)-ra.steady_state(j));
				unit = 'Percentage points';
			else
				x = 100*(oo.endo_simul(i,:)/oo.steady_state(i)-1);
				y = 100*(ra.endo_simul(j,:)/ra.steady_state(j)-1);
				unit = '% from steady state';
			end
			plot(time,x,'LineWidth',2);
			hold on;
			plot(time,y,'--','LineWidth',2);
			ylabel(unit);
			grid on;
			xlim([0 40]);
			xlabel('Years');
			title(titles{q});
			if q == 1
				legend('OLG','RANK','Location','southeast');
			end
	end
end
