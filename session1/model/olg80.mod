// Pedagogical 80-cohort OLG model for Session 1.
// Ages 20--99, flexible prices, no nominal rigidities.
// The equilibrium interest rate closes the capital/asset market.

@#define J = 80
@#ifndef SHOCK_TYPE
  @#define SHOCK_TYPE = 1 // 1: TFP; 2: depreciation (financial crisis)
@#endif

var Z dep y w r k C
@#for j in 1:J
    c_@{j}
@#endfor
@#for j in 1:J-1
    a_@{j}
@#endfor
;

varexo epsZ epsD;

parameters beta sigma alpha delta rhoZ rhoD Lbar
@#for j in 1:J
    z_@{j} N_@{j}
@#endfor
;

// Calibration (identical in olg80.mod and rank80.mod; olg80_bc.mod adds a transfer process).
beta  = 0.97;
sigma = 2.00;
alpha = 0.33;
delta = 0.06;
rhoZ  = 0.90;
rhoD  = 0.90;

// Exogenous efficiency profile: hump-shaped earnings until age 64,
// retirement thereafter. Cohort masses are equal in this first version.
@#for j in 1:J
    N_@{j} = 1/@{J};
    @#if j <= 45
        z_@{j} = 0.35 + 0.85*exp(-0.5*((@{j+19}-50)/14)^2);
    @#else
        z_@{j} = 0;
    @#endif
@#endfor

Lbar = 0
@#for j in 1:J
    + N_@{j}*z_@{j}
@#endfor
;

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

steady(maxit=1000, solve_algo=4);
check;

@#if SHOCK_TYPE == 1
    shocks;
        var epsZ;
        periods 1;
        values -0.01;
    end;
@#else
    shocks;
        var epsD;
        periods 1;
        values 0.02;
    end;
@#endif

perfect_foresight_setup(periods=100);
perfect_foresight_solver(maxit=200, tolf=1e-9);

assert(oo_.deterministic_simulation.status == 1, 'Perfect-foresight solution failed');
if ~isfolder('results'), mkdir('results'); end

@#if SHOCK_TYPE == 1
    save('results/results_tfp.mat','M_','oo_','options_');
@#else
    save('results/results_dep.mat','M_','oo_','options_');
@#endif

// ---------------------------------------------------------------------
// Figures of the handout: calibration and stationary equilibrium.
// Transition figures (OLG against RANK) are drawn at the end of rank80.mod.
// ---------------------------------------------------------------------
@#if SHOCK_TYPE == 1
verbatim;
figdir = fullfile('..','figures'); if ~isfolder(figdir), mkdir(figdir); end
names = cellstr(M_.endo_names); pn = cellstr(M_.param_names);
J = @{J}; age = (20:19+J)';
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
Lb = M_.params(strcmp(pn,'Lbar')); al = M_.params(strcmp(pn,'alpha'));
be = M_.params(strcmp(pn,'beta')); de = M_.params(strcmp(pn,'delta')); si = M_.params(strcmp(pn,'sigma'));
rr = linspace(-0.015,0.04,120); A = nan(size(rr)); Kd = A; h = (0:J-1)';
for i = 1:numel(rr)
    R = 1+rr(i); Kd(i) = (al*Lb^(1-al)/(rr(i)+de))^(1/(1-al)); w = (1-al)*Kd(i)^al*Lb^(-al);
    g = (be*R)^(1/si); c = sum(w*z./R.^h)/sum(g.^h./R.^h)*g.^h;
    a = zeros(J,1); a(1) = w*z(1)-c(1);
    for j = 2:J, a(j) = R*a(j-1)+w*z(j)-c(j); end
    A(i) = Nj'*a;
end
figure('Color','w','Position',[100 100 620 430]); hold('on');
plot(A,100*rr,'Color',blue,'LineWidth',2.2); plot(Kd,100*rr,'k','LineWidth',1.6);
plot(ss('k'),100*ss('r'),'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',7);
grid('on'); xlabel('Assets / capital per capita'); ylabel('Net interest rate r (%)'); xlim([0 12]);
legend('Asset supply','Capital demand','Equilibrium','Location','southwest');
exportgraphics(gcf,fullfile(figdir,'olg80_asset_market.pdf'),'ContentType','vector');
end;
@#endif
