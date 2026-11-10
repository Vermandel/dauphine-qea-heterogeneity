// Pedagogical 80-cohort OLG model for Session 1.
// Ages 20--99, flexible prices, no nominal rigidities.
// The equilibrium interest rate closes the capital/asset market.

@#define J = 80
@#ifndef SHOCK_TYPE
  @#define SHOCK_TYPE = 1 // 1: TFP; 2: depreciation (financial crisis)
@#endif

var
	Z		// Aggregate total factor productivity.
	dep		// Time-varying depreciation rate.
	y		// Aggregate output.
	w		// Wage per efficiency unit.
	r		// Net return on capital.
	k		// End-of-period aggregate capital stock.
	C		// Aggregate consumption.
@#for j in 1:J
	c_@{j}		// Consumption of age group @{j} (age @{j+19}).
@#endfor
@#for j in 1:J-1
	a_@{j}		// End-of-period assets of age group @{j}.
@#endfor
;

varexo
	epsZ		// Innovation to aggregate productivity.
	epsD		// Innovation to the depreciation rate.
;

parameters
	beta		// Annual subjective discount factor.
	sigma		// CRRA coefficient (inverse intertemporal elasticity).
	alpha		// Capital share in output.
	delta		// Baseline annual depreciation rate.
	rhoZ		// Persistence of aggregate productivity.
	rhoD		// Persistence of depreciation.
	Lbar		// Aggregate labor in efficiency units.
@#for j in 1:J
	z_@{j}		// Labor efficiency of age group @{j}.
	N_@{j}		// Population share of age group @{j}.
@#endfor
;

// Calibration (identical in olg80.mod and rank80.mod; olg80_bc.mod adds a transfer process).
beta  = 0.97;	// Annual subjective discount factor.
sigma = 2.00;	// CRRA coefficient (inverse intertemporal elasticity).
alpha = 0.33;	// Capital share in production.
delta = 0.06;	// Baseline annual depreciation rate.
rhoZ  = 0.90;	// Persistence of aggregate productivity.
rhoD  = 0.90;	// Persistence of depreciation shocks.

// Cohort j is age j+19. Working efficiency is a base of 0.35 plus a smooth
// mid-career bonus of up to 0.85 at age 50; 14 years controls its spread.
// At 65 retirement sets efficiency to zero. Equal cohort masses reflect
// certain survival until age 99 and zero population growth.
@#for j in 1:J
	N_@{j} = 1/@{J};	// Equal population mass at every age.
	@#if j <= 45
		z_@{j} = 0.35 + 0.85*exp(-0.5*((@{j+19}-50)/14)^2);	// Base plus age-50 bonus.
	@#else
		z_@{j} = 0;	// No labor efficiency after retirement.
	@#endif
@#endfor

// Each age group supplies its labor efficiency times its population share.
// Adding these contributions gives the total labor input used by firms.
Lbar = 0
@#for j in 1:J
	+ N_@{j}*z_@{j}
@#endfor
;	// Labor supply weighted by cohort mass and efficiency.

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
    c_1 + a_1 = w*z_1;

    @#for j in 2:J-1
        [name='Budget cohort @{j}']
        c_@{j} + a_@{j} = (1+r)*a_@{j-1}(-1) + w*z_@{j};
    @#endfor

    [name='Terminal budget age 99']
    c_@{J} = (1+r)*a_@{J-1}(-1) + w*z_@{J};

    @#for j in 1:J-1
        [name='Euler cohort @{j}']
        c_@{j}^(-sigma) = beta*(1+r(+1))*c_@{j+1}(+1)^(-sigma);
    @#endfor
end;

// The household life cycle is analytic at a candidate interest rate.
// Only the rate that clears the aggregate asset market requires a scalar root.
steady_state_model;
    Z   = 1;
    dep = delta;
    r   = olg80_interest_rate(alpha,beta,delta,sigma,Lbar);
    k   = (alpha*Lbar^(1-alpha)/(r+delta))^(1/(1-alpha));
    y   = k^alpha*Lbar^(1-alpha);
    w   = (1-alpha)*y/Lbar;

    R = 1+r;                         // Gross return.
    g = (beta*R)^(1/sigma);          // Euler growth: c_j = g^(j-1)*c_1.
    // Lifetime resources divided by discounted consumption units determine c_1.
    c_1 = w*(0
    @#for j in 1:J
        + z_@{j}/R^@{j-1}
    @#endfor
    )/(0
    @#for j in 1:J
        + (g/R)^@{j-1}
    @#endfor
    );
    @#for j in 2:J
        c_@{j} = g*c_@{j-1};
    @#endfor

    a_1 = w*z_1-c_1;               // Youngest cohort starts with no assets.
    @#for j in 2:J-1
        a_@{j} = R*a_@{j-1}+w*z_@{j}-c_@{j};
    @#endfor
    C = 0
    @#for j in 1:J
        + N_@{j}*c_@{j}
    @#endfor
    ;
end;

steady(maxit=1000, solve_algo=4);
check;

@#if SHOCK_TYPE == 1
    shocks;
		var epsZ;	// Shock to aggregate productivity.
		periods 1;	// Occurs in the first simulated period.
		values -0.01;	// Log productivity falls by 0.01 (about 1%).
    end;
@#else
    shocks;
		var epsD;	// Shock to the depreciation rate.
		periods 1;	// Occurs in the first simulated period.
		values 0.02;	// Depreciation rises by 2 percentage points.
    end;
@#endif

perfect_foresight_setup(periods=100);
perfect_foresight_solver(maxit=200, tolf=1e-8);
// Dynare prepares the chosen shock and an initial guess over 100 periods.
// It then solves for consumption, saving, capital and prices along that path,
// with households knowing the future shock sequence when they choose.

assert(oo_.deterministic_simulation.status == 1, 'Perfect-foresight solution failed');

@#if SHOCK_TYPE == 1
    save('results_tfp.mat','M_','oo_','options_');
@#else
    save('results_dep.mat','M_','oo_','options_');
@#endif

// ---------------------------------------------------------------------
// Calibration and stationary equilibrium (subplots).
// The OLG-vs-RANK transition figures are drawn at the end of rank80.mod.
// ---------------------------------------------------------------------
@#if SHOCK_TYPE == 1
verbatim;
% Prepare age-indexed stationary profiles and the teaching figures.
names = cellstr(M_.endo_names);
pn    = cellstr(M_.param_names);
J     = @{J};
age   = (20:19+J)';

% Steady-state objects by age
ss  = @(n) oo_.steady_state(strcmp(names,n));
z   = arrayfun(@(j) M_.params(strcmp(pn,sprintf('z_%d',j))), 1:J)';
css = arrayfun(@(j) ss(sprintf('c_%d',j)), 1:J)';
ass = [arrayfun(@(j) ss(sprintf('a_%d',j)), 1:J-1)'; 0];

% No borrowing limit: distribute one unit of extra cash over remaining life.
% The household path is analytic at stationary prices; only the market-clearing rate needs a scalar root.
R    = 1+ss('r');	% Gross stationary return.
beta = M_.params(strcmp(pn,'beta'));	% Annual household discount factor.
sigma = M_.params(strcmp(pn,'sigma'));	% Relative risk aversion.
g    = (beta*R)^(1/sigma);	% Euler consumption growth factor.
mpc  = arrayfun(@(j) 1/sum((g/R).^(0:J-j)), (1:J)');	% Marginal consumption from today's windfall.
a_ss = ass;	% End-of-period assets in the unconstrained stationary equilibrium.
save('mpc_bc0.mat','age','mpc','a_ss');

figure('Name','OLG steady state');

subplot(2,2,1);
	plot(age,z,'LineWidth',2);
	grid on;
	xlabel('Age');
	title('Labor efficiency e_j');

subplot(2,2,2);
	plot(age,ss('w')*z,'LineWidth',2);
	hold on;
	plot(age,css,'--','LineWidth',2);
	grid on;
	xlabel('Age');
	title('Labor income and consumption');
	legend('Labor income','Consumption','Location','east');

subplot(2,2,3);
	plot(age,css,'LineWidth',2);
	grid on;
	xlabel('Age');
	title('Consumption profile');

subplot(2,2,4);
	plot(age,ass,'LineWidth',2);
	hold on;
	yline(0,'k:');
	grid on;
	xlabel('Age');
	title('Wealth profile');
end;
@#endif

// Responses of the OLG aggregates after the selected productivity or depreciation shock.
verbatim;
	% Compare each aggregate with its own stationary value over the first 40 years.
	names  = cellstr(M_.endo_names);	% Names of endogenous variables.
	idx    = @(n) find(strcmp(names,n),1);	% Position of an aggregate.
	time   = 0:size(oo_.endo_simul,2)-1;	% Years since the shock.
	vars   = {'C','y','k','r'};	% Aggregates to display.
	labels = {'Consumption','Output','Capital','Net interest rate'};	% Panel titles.

@#if SHOCK_TYPE == 1
	figure('Name','OLG: productivity shock');
@#else
	figure('Name','OLG: depreciation shock');
@#endif
	for q = 1:numel(vars)
		subplot(2,2,q);
			i = idx(vars{q});	% Dynare row of the selected aggregate.
			if strcmp(vars{q},'r')
				response = 100*(oo_.endo_simul(i,:)-oo_.steady_state(i));
				unit = 'Percentage points';
			else
				response = 100*(oo_.endo_simul(i,:)/oo_.steady_state(i)-1);
				unit = '% from steady state';
			end
			plot(time,response,'LineWidth',2);
			ylabel(unit);
			yline(0,'k:');
			grid on;
			xlim([0 40]);
			xlabel('Years');
			title(labels{q});
	end
end;
