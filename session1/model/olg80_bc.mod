// Appendix model: olg80.mod plus the borrowing limit a_j >= 0 (BC=1) and a
// balanced transfer experiment (SHOCK_TYPE=2, windfall to cohort TARGET).
// Options: -DBC=0/1  -DSHOCK_TYPE=1/2  -DTARGET=<cohort>  -DHORIZON=<periods>

@#define J = 80
@#ifndef BC
  @#define BC = 1
@#endif
@#ifndef HORIZON
  @#define HORIZON = 200
@#endif
@#ifndef SHOCK_TYPE
  @#define SHOCK_TYPE = 1 // 1: TFP; 2: balanced redistribution
@#endif
@#ifndef TARGET
  @#define TARGET = 31   // cohort index; age = TARGET + 19
@#endif

var Z T y w r k C
@#for j in 1:J
    c_@{j}
@#endfor
@#for j in 1:J-1
    a_@{j}
@#endfor
;

varexo epsZ epsT;

parameters beta sigma alpha delta rhoZ rhoT Lbar constrained
@#for j in 1:J
    z_@{j} N_@{j} psi_@{j}
@#endfor
;

// Calibration (as in olg80.mod, with a transfer process rhoT instead of the depreciation process rhoD).
beta  = 0.97;
sigma = 2.00;
alpha = 0.33;
delta = 0.06;
rhoZ  = 0.90;
rhoT  = 0.00;

// Exogenous efficiency profile: hump-shaped earnings until age 64,
// retirement thereafter. Cohort masses are equal in this first version.
@#for j in 1:J
    N_@{j} = 1/@{J};
    @#if j <= 45
        z_@{j} = 0.35 + 0.85*exp(-0.5*((@{j+19}-50)/14)^2);
    @#else
        z_@{j} = 0;
    @#endif
    psi_@{j} = -1/@{J};
    @#if j == TARGET
        psi_@{j} = 1-1/@{J};
    @#endif
@#endfor

Lbar = 0
@#for j in 1:J
    + N_@{j}*z_@{j}
@#endfor
;
constrained = @{BC};

model;
    [name='TFP law of motion']
    log(Z) = rhoZ*log(Z(-1)) + epsZ;

    [name='Balanced redistribution shock']
    T = rhoT*T(-1) + epsT;

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
    c_1 + a_1 = w*z_1 + psi_1*T;

    @#for j in 2:J-1
        [name='Budget cohort @{j}']
        c_@{j} + a_@{j} = (1+r)*a_@{j-1}(-1) + w*z_@{j} + psi_@{j}*T;
    @#endfor

    [name='Terminal budget age 99']
    c_@{J} = (1+r)*a_@{J-1}(-1) + w*z_@{J} + psi_@{J}*T;

    @#for j in 1:J-1
        [name='Euler cohort @{j}']
        @#if BC == 1
        // ASV borrowing limit, with no portfolio adjustment costs:
        // c <= cash on hand is equivalent to end-of-period a >= 0.
        c_@{j} = min(c_@{j+1}(+1)/(beta*(1+r(+1)))^(1/sigma),
            w*z_@{j}+psi_@{j}*T
            @#if j > 1
              +(1+r)*a_@{j-1}(-1)
            @#endif
        );
        @#else
        c_@{j}^(-sigma) = beta*(1+r(+1))*c_@{j+1}(+1)^(-sigma);
        @#endif
    @#endfor
end;

steady(maxit=1000, solve_algo=4);

@#if SHOCK_TYPE == 1
    shocks;
        var epsZ;
        periods 1;
        values -0.01;
    end;
@#else
    shocks;
        var epsT;
        periods 1;
        values 0.0001;
    end;
@#endif

perfect_foresight_setup(periods=@{HORIZON});
perfect_foresight_solver(maxit=200, tolf=1e-10);

assert(oo_.deterministic_simulation.status == 1, 'Perfect-foresight solution failed');
if ~isfolder('results'), mkdir('results'); end

@#if SHOCK_TYPE == 1
    save('results/results_bc@{BC}_tfp.mat','M_','oo_','options_');
@#else
    save('results/results_bc@{BC}_age@{TARGET+19}.mat','M_','oo_','options_');
@#endif

// ---------------------------------------------------------------------
// Output. SHOCK_TYPE=2: general-equilibrium response to the windfall.
// SHOCK_TYPE=1: TFP figure (needs olg80 and rank80 runs) and MPC figure.
// ---------------------------------------------------------------------
@#if SHOCK_TYPE == 2
verbatim;
names = cellstr(M_.endo_names); pn = cellstr(M_.param_names);
ic = find(strcmp(names,'c_@{TARGET}'),1);
psi = M_.params(strcmp(pn,'psi_@{TARGET}'));
fprintf('Borrowing limit %d, age %d: consumption response / windfall = %.6f\n', ...
    @{BC}, @{TARGET}+19, (oo_.endo_simul(ic,2)-oo_.steady_state(ic))/(psi*1e-4));
end;
@#else
verbatim;
figdir = fullfile('..','figures'); if ~isfolder(figdir), mkdir(figdir); end
names = cellstr(M_.endo_names); pn = cellstr(M_.param_names);
J = @{J}; age = (20:19+J)'; t = 0:size(oo_.endo_simul,2)-1;
blue = [0 0.30 0.60]; orange = [0.85 0.45 0.05]; red = [0.80 0.15 0.10]; grey = [0.45 0.45 0.45 0.22];
z = arrayfun(@(j) M_.params(strcmp(pn,sprintf('z_%d',j))), 1:J)';
pdev = @(o,i) 100*(o.endo_simul(i,:)/o.steady_state(i)-1);

% MPC out of a one-off windfall at each age, prices held fixed.
step = 1e-4; mpc = zeros(J,2); aprof = zeros(J,2);
for bc = 0:1
    M2 = M_; M2.params(strcmp(pn,'constrained')) = bc;
    ys = olg80_bc_steadystate([],[],M2,[]);
    s.r = ys(strcmp(names,'r')); s.w = ys(strcmp(names,'w'));
    s.c = arrayfun(@(j) ys(strcmp(names,sprintf('c_%d',j))), 1:J)';
    s.a = [arrayfun(@(j) ys(strcmp(names,sprintf('a_%d',j))), 1:J-1)'; 0];
    S(bc+1) = s; aprof(:,bc+1) = s.a;
    be = M_.params(strcmp(pn,'beta')); si = M_.params(strcmp(pn,'sigma'));
    prev = [0; s.a(1:end-1)];
    for j = 1:J
        inc = s.w*z(j:end); c0 = lifecycle_bc(inc,1+s.r,be,si,prev(j),bc);
        inc(1) = inc(1)+step; c1 = lifecycle_bc(inc,1+s.r,be,si,prev(j),bc);
        mpc(j,bc+1) = (c1(1)-c0(1))/step;
    end
end
[cf,af] = lifecycle_bc(S(1).w*z,1+S(1).r,be,si,0,1); prev = [0; af(1:end-1)]; fixed = zeros(J,1);
for j = 1:J
    inc = S(1).w*z(j:end); c0 = lifecycle_bc(inc,1+S(1).r,be,si,prev(j),1);
    inc(1) = inc(1)+step; c1 = lifecycle_bc(inc,1+S(1).r,be,si,prev(j),1);
    fixed(j) = (c1(1)-c0(1))/step;
end
figure('Color','w','Position',[50 50 1100 420]); tiledlayout(1,2,'TileSpacing','compact');
nexttile; plot(age,aprof(:,1),'LineWidth',1.8); hold('on'); plot(age,aprof(:,2),'LineWidth',1.8);
yline(0,':k'); grid('on'); xlabel('Age'); ylabel('End-of-period assets'); title('Stationary asset profiles');
legend('No borrowing limit','ASV limit: a >= 0','Location','northwest');
nexttile; plot(age,mpc(:,1),'LineWidth',1.8); hold('on'); plot(age,mpc(:,2),'LineWidth',1.8);
plot(age,fixed,'--','LineWidth',1.2); grid('on'); ylim([0 1.05]); xlabel('Age'); ylabel('Impact MPC');
title('One-off windfall, prices held fixed');
legend('No limit: own steady state','ASV limit: own steady state','ASV limit: baseline prices','Location','north');
exportgraphics(gcf,fullfile(figdir,'olg80_borrowing_mpc.pdf'),'ContentType','vector');
fprintf('Binding ages: %s\n',mat2str(age(aprof(1:J,2)<1e-8 & age<age(end))'));
fprintf('MPC at ages 25, 35, 75, 99 (no limit): %s\n',mat2str(mpc([6 16 56 80],1)',6));
fprintf('MPC at ages 25, 35, 75, 99 (limit):    %s\n',mat2str(mpc([6 16 56 80],2)',6));

% TFP transition with and without the limit, against RANK.
if @{BC} == 1
    O = load(fullfile('results','results_tfp.mat')); R = load(fullfile('results','results_rank_tfp.mat'));
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
end
end;
@#endif
