// Borrowing-constrained OLG: a_j >= 0 at each age.
// For the unconstrained analytic life cycle, run olg80.mod instead.
// Option: -DHORIZON=<periods>

@#define J = 80
@#ifndef HORIZON
  @#define HORIZON = 200
@#endif

var
	Z		// Aggregate productivity.
	y		// Aggregate output.
	w		// Wage per efficiency unit.
	r		// Net capital return.
	k		// End-of-period capital stock.
	C		// Aggregate consumption.
@#for j in 1:J
	c_@{j}		// Consumption of age group @{j}.
@#endfor
@#for j in 1:J-1
	a_@{j}		// End-of-period assets of age group @{j}.
@#endfor
;

varexo
	epsZ		// Aggregate productivity innovation.
	epsW		// Household windfall in the MPC experiment.
;

parameters
	beta		// Annual household discount factor.
	sigma		// Relative risk aversion (inverse intertemporal elasticity).
	alpha		// Capital share of output.
	delta		// Annual depreciation rate.
	rhoZ		// Persistence of productivity.
	Lbar		// Aggregate efficient labor.
@#for j in 1:J
	z_@{j}		// Labor efficiency of age group @{j}.
	N_@{j}		// Population share of age group @{j}.
	dW_@{j}		// Selects the age group receiving a windfall.
@#endfor
;

// Calibration (as in olg80.mod).
beta  = 0.97;	// Annual subjective discount factor.
sigma = 2.00;	// Relative risk aversion (inverse intertemporal elasticity).
alpha = 0.33;	// Capital share in production.
delta = 0.06;	// Annual depreciation rate.
rhoZ  = 0.90;	// Persistence of productivity shocks.

// Cohort j is age j+19. Working efficiency is 0.35 plus a smooth bonus
// of up to 0.85 at age 50; retirement sets it to zero from age 65.
// Equal cohort masses imply the same population weight at each age.
@#for j in 1:J
	N_@{j} = 1/@{J};	// Equal population mass at every age.
	dW_@{j} = 0;	// No cohort receives the windfall initially.
	@#if j <= 45
		z_@{j} = 0.35 + 0.85*exp(-0.5*((@{j+19}-50)/14)^2);	// Base plus age-50 bonus.
	@#else
		z_@{j} = 0;	// No labor supply after retirement.
	@#endif
@#endfor

// Each age group contributes its efficiency times its population share.
// Adding them gives the total labor input used by firms.
Lbar = 0
@#for j in 1:J
	+ N_@{j}*z_@{j}
@#endfor
;	// Aggregate efficient labor.

model;
    [name='TFP law of motion']
    log(Z) = rhoZ*log(Z(-1)) + epsZ;

    [name='Production']
    y = Z*k(-1)^alpha*Lbar^(1-alpha);

    [name='Wage']
    w = (1-alpha)*y/Lbar;

    [name='Net return on capital']
    r = alpha*y/k(-1)-delta;

    [name='Aggregate consumption']
    C = 0
    @#for j in 1:J
        + N_@{j}*c_@{j}
    @#endfor
    ;

    [name='Asset-market clearing']
    k = 0
    @#for j in 1:J-1
        + N_@{j}*a_@{j}
    @#endfor
    ;

    [name='Budget age 20']
    c_1 + a_1 = w*z_1 + dW_1*epsW;

    @#for j in 2:J-1
        [name='Budget cohort @{j}']
        c_@{j} + a_@{j} = (1+r)*a_@{j-1}(-1) + w*z_@{j} + dW_@{j}*epsW;
    @#endfor

    [name='Terminal budget age 99']
    c_@{J} = (1+r)*a_@{J-1}(-1) + w*z_@{J} + dW_@{J}*epsW;

    @#for j in 1:J-1
        [name='Euler cohort @{j}']
        // Cash on hand caps consumption, ensuring end-of-period assets a_j >= 0.
        c_@{j} = min(c_@{j+1}(+1)/(beta*(1+r(+1)))^(1/sigma),
            w*z_@{j} + dW_@{j}*epsW
            @#if j > 1
              +(1+r)*a_@{j-1}(-1)
            @#endif
        );
    @#endfor
end;

steady(maxit=1000, solve_algo=4);

shocks;
	var epsZ;	// Productivity innovation.
	periods 1;	// First simulated period.
	values -0.01;	// Log productivity falls by about 1%.
end;

perfect_foresight_setup(periods=@{HORIZON});
perfect_foresight_solver(maxit=200, tolf=1e-8);
// Prepare the productivity shock over the chosen horizon and solve the
// consumption, saving, capital and price paths under perfect foresight.

assert(oo_.deterministic_simulation.status == 1, 'Perfect-foresight solution failed');

save('results_bc1_tfp.mat','M_','oo_','options_');

// ---------------------------------------------------------------------
// Fixed-price MPC: an unexpected windfall changes one cohort's current cash.
// Keep inherited assets, wages and interest rates at their own steady state.
// The dW_j / epsW terms remain available for separate GE experiments.
// ---------------------------------------------------------------------
verbatim;
% Compute constrained MPCs from remaining-life prefix budgets at fixed prices.
names = cellstr(M_.endo_names);
pn    = cellstr(M_.param_names);
J     = @{J};
step  = 1e-3;
ss    = @(n) oo_.steady_state(strcmp(names,n));
a_ss  = [arrayfun(@(j) ss(sprintf('a_%d',j)), 1:J-1)'; 0];
c_ss  = arrayfun(@(j) ss(sprintf('c_%d',j)), 1:J)';
z     = arrayfun(@(j) M_.params(strcmp(pn,sprintf('z_%d',j))), 1:J)';
w_ss  = ss('w');
R_ss  = 1+ss('r');
beta  = M_.params(strcmp(pn,'beta'));
sigma = M_.params(strcmp(pn,'sigma'));
mpc   = zeros(J,1);
for j = 1:J
	inherited = 0;
	if j > 1
		inherited = a_ss(j-1);
	end
	income = w_ss*z(j:J);
	income(1) = income(1)+step;
	[c_windfall,~] = lifecycle_bc(income,R_ss,beta,sigma,inherited);
	mpc(j) = (c_windfall(1)-c_ss(j))/step;
end
end;

// ---------------------------------------------------------------------
// Output: stationary assets and MPC profiles (subplots).
// ---------------------------------------------------------------------
verbatim;
age = (20:19+@{J})';
save('mpc_bc1.mat','age','mpc','a_ss');

figure('Name','Borrowing limit');

subplot(1,2,1);
	plot(age,a_ss,'LineWidth',1.8);
	yline(0,':k');
	grid on;
	xlabel('Age');
	title('Stationary asset profile');

subplot(1,2,2);
	plot(age,mpc,'LineWidth',1.8);
	grid on;
	ylim([0 1.05]);
	xlabel('Age');
	title('Impact MPC out of a windfall');

fprintf('MPC at ages 25, 35, 75, 99 (borrowing limit): %s\n',mat2str(mpc([6 16 56 80])',6));
end;
