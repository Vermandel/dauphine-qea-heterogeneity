// Representative-agent (RANK/RBC) benchmark for the 80-cohort OLG.
// Same technology, labor input and TFP shock as olg80.mod. The discount
// factor is set so that the RA steady state reproduces the OLG interest
// rate: both economies then share K*, Y*, w* and C*, and differ only in
// how households respond to the shock.
@#define J = 80
@#ifndef SHOCK_TYPE
  @#define SHOCK_TYPE = 1 // 1: TFP; 2: depreciation (financial crisis)
@#endif
var
	Z		// Aggregate productivity.
	dep		// Time-varying depreciation rate.
	y		// Aggregate output.
	w		// Wage per efficiency unit.
	r		// Net capital return.
	k		// End-of-period capital stock.
	C		// Aggregate consumption.
;
varexo
	epsZ		// Aggregate productivity innovation.
	epsD		// Depreciation-rate innovation.
;
parameters
	beta		// Annual household discount factor in the OLG calibration.
	sigma		// Relative risk aversion (inverse intertemporal elasticity).
	alpha		// Capital share of output.
	delta		// Baseline annual depreciation rate.
	rhoZ		// Persistence of productivity.
	rhoD		// Persistence of depreciation.
	Lbar		// Aggregate efficient labor supplied.
	betaRA		// RA discount factor chosen to match the OLG return.
@#for j in 1:J
	z_@{j}		// Labor efficiency of age group @{j}.
	N_@{j}		// Population share of age group @{j}.
@#endfor
;
// Calibration (identical in olg80.mod and rank80.mod; olg80_bc.mod adds a transfer process).
beta  = 0.97;	// OLG annual discount factor.
sigma = 2.00;	// Relative risk aversion (inverse intertemporal elasticity).
alpha = 0.33;	// Capital share in production.
delta = 0.06;	// Baseline annual depreciation rate.
rhoZ  = 0.90;	// Persistence of productivity shocks.
rhoD  = 0.90;	// Persistence of depreciation shocks.

// Cohort j is age j+19. Working efficiency is a base of 0.35 plus a smooth
// age-50 bonus of up to 0.85; 14 years controls its spread. From age 65,
// retirement sets efficiency to zero. Cohort masses are equal.
@#for j in 1:J
	N_@{j} = 1/@{J};	// Equal population mass at each age.
	@#if j <= 45
		z_@{j} = 0.35 + 0.85*exp(-0.5*((@{j+19}-50)/14)^2);	// Base plus age-50 bonus.
	@#else
		z_@{j} = 0;	// No labor supply after retirement.
	@#endif
@#endfor

// Add each age group's efficiency weighted by its population share.
// The sum is the total labor input used in production.
Lbar = 0
@#for j in 1:J
	+ N_@{j}*z_@{j}
@#endfor
;	// Aggregate efficient labor.
// Match the representative household's discount factor to the saved OLG return.
@#if SHOCK_TYPE == 1
olg_ = load('results_tfp.mat','oo_','M_');	// OLG economy after the productivity shock.
@#else
olg_ = load('results_dep.mat','oo_','M_');	// OLG economy after the depreciation shock.
@#endif
rOLG = olg_.oo_.steady_state(strcmp(cellstr(olg_.M_.endo_names),'r'));	// OLG interest rate.
set_param_value('betaRA',1/(1+rOLG));	// Same stationary return in the RA economy.

model;
    [name='TFP law of motion']
    log(Z) = rhoZ*log(Z(-1)) + epsZ;
    [name='Depreciation process']
    dep = (1-rhoD)*delta + rhoD*dep(-1) + epsD;

    [name='Production']
    y = Z*k(-1)^alpha*Lbar^(1-alpha);
    [name='Wage']
    w = (1-alpha)*y/Lbar;
    [name='Net return on capital']
    r = alpha*y/k(-1)-dep;
    [name='Representative budget']
    C + k = (1+r)*k(-1) + w*Lbar;
    [name='Representative Euler equation']
    C^(-sigma) = betaRA*(1+r(+1))*C(+1)^(-sigma);
end;

steady_state_model;
    Z = 1;
    dep = delta;
    r = 1/betaRA - 1;
    k = (alpha*Lbar^(1-alpha)/(r+delta))^(1/(1-alpha));
    y = k^alpha*Lbar^(1-alpha);
    w = (1-alpha)*y/Lbar;
    C = y - delta*k;
end;
steady;
check;

@#if SHOCK_TYPE == 1
    shocks;
		var epsZ;	// Shock to productivity.
		periods 1;	// First simulated period.
		values -0.01;	// Log productivity falls by about 1%.
    end;
@#else
    shocks;
		var epsD;	// Shock to depreciation.
		periods 1;	// First simulated period.
		values 0.02;	// Depreciation rises by 2 percentage points.
    end;
@#endif
perfect_foresight_setup(periods=100);
perfect_foresight_solver(maxit=200, tolf=1e-8);
// Prepare the same 100-period shock path as OLG, then solve all RA choices
// and market prices along it under perfect foresight.
assert(oo_.deterministic_simulation.status == 1, 'RANK perfect-foresight solution failed');
@#if SHOCK_TYPE == 1
    save('results_rank_tfp.mat','M_','oo_','options_');
@#else
    save('results_rank_dep.mat','M_','oo_','options_');
@#endif

// ---------------------------------------------------------------------
// Transition figure (subplots): OLG (saved by olg80.mod) against RANK, same shock.
// ---------------------------------------------------------------------
verbatim;
% Compare the saved OLG transition with the current RA path.
@#if SHOCK_TYPE == 1
	O = load('results_tfp.mat');
@#else
	O = load('results_dep.mat');
@#endif
oo = O.oo_;		% OLG simulated transition.
on = cellstr(O.M_.endo_names);	% OLG variable names.
rn = cellstr(M_.endo_names);	% RA variable names.
io = @(n) find(strcmp(on,n),1);	% Position of an OLG variable.
ir = @(n) find(strcmp(rn,n),1);	% Position of an RA variable.
t  = 0:size(oo.endo_simul,2)-1;	% Simulation periods.

% Percentage deviation from steady state (percentage points for r)
pdev = @(o,i) 100*(o.endo_simul(i,:)/o.steady_state(i)-1);
ppt  = @(o,i) 100*(o.endo_simul(i,:)-o.steady_state(i));

vars = {'C','y','k','r'};
ttl  = {'Consumption','Output','Capital','Net interest rate'};

figure('Name','OLG vs RANK');
for q = 1:4
	subplot(2,2,q);
	hold on;
	if q == 4
		f = ppt;
		unit = 'Percentage points';
	else
		f = pdev;
		unit = '% from steady state';
	end
	plot(t,f(oo,io(vars{q})),'LineWidth',2);
	plot(t,f(oo_,ir(vars{q})),'--','LineWidth',2);
	ylabel(unit);
	title(ttl{q});
	xlabel('Years');
	grid on;
	xlim([0 40]);
	if q == 1
		legend('OLG','RANK','Location','southeast');
	end
end
end;
