// Representative-agent (RANK/RBC) benchmark for the 80-cohort OLG.
// Same technology, labor input and TFP shock as olg80.mod. The discount
// factor is set so that the RA steady state reproduces the OLG interest
// rate: both economies then share K*, Y*, w* and C*, and differ only in
// how households respond to the shock.
@#define J = 80
@#ifndef SHOCK_TYPE
  @#define SHOCK_TYPE = 1 // 1: TFP; 2: depreciation (financial crisis)
@#endif
var Z dep y w r k C;
varexo epsZ epsD;
parameters beta sigma alpha delta rhoZ rhoD Lbar betaRA
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
olg_ = load(fullfile('results','results_tfp.mat'),'oo_','M_');
rOLG = olg_.oo_.steady_state(strcmp(cellstr(olg_.M_.endo_names),'r'));
set_param_value('betaRA',1/(1+rOLG));

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
assert(oo_.deterministic_simulation.status == 1, 'RANK perfect-foresight solution failed');
@#if SHOCK_TYPE == 1
    save('results/results_rank_tfp.mat','M_','oo_','options_');
@#else
    save('results/results_rank_dep.mat','M_','oo_','options_');
@#endif

// ---------------------------------------------------------------------
// Transition figure: OLG (saved by olg80.mod) against RANK, same shock.
// ---------------------------------------------------------------------
verbatim;
figdir = fullfile('..','figures'); if ~isfolder(figdir), mkdir(figdir); end
@#if SHOCK_TYPE == 1
    O = load(fullfile('results','results_tfp.mat')); figname = 'olg80_tfp_transition.pdf';
@#else
    O = load(fullfile('results','results_dep.mat')); figname = 'olg80_depreciation.pdf';
@#endif
oo = O.oo_; J = @{J}; age = (20:19+J)';
on = cellstr(O.M_.endo_names); rn = cellstr(M_.endo_names);
io = @(n) find(strcmp(on,n),1); ir = @(n) find(strcmp(rn,n),1);
pdev = @(o,i) 100*(o.endo_simul(i,:)/o.steady_state(i)-1);
ppt  = @(o,i) 100*(o.endo_simul(i,:)-o.steady_state(i));
t = 0:size(oo.endo_simul,2)-1;
grey = [0.45 0.45 0.45 0.22]; blue = [0 0.30 0.60]; red = [0.80 0.15 0.10];

figure('Color','w','Position',[100 100 1000 600]);
@#if SHOCK_TYPE == 1
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile([2 1]);
@#else
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
@#endif
hold('on');
for j = 1:J, plot(t,pdev(oo,io(sprintf('c_%d',j))),'Color',grey,'LineWidth',0.8); end
h1 = plot(t,pdev(oo,io('C')),'Color',blue,'LineWidth',2.5);
h2 = plot(t,pdev(oo_,ir('C')),'--','Color',red,'LineWidth',2);
h0 = plot(nan,nan,'Color',grey(1:3),'LineWidth',0.8);
title('Consumption'); ylabel('% deviation'); grid('on'); xlim([0 40]);
legend([h0 h1 h2],{'OLG by age','OLG aggregate','RANK'},'Location','southeast');
@#if SHOCK_TYPE == 1
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
@#else
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
@#endif
xlabel(tl,'Years (unexpected shock at date 1)');
exportgraphics(gcf,fullfile(figdir,figname),'ContentType','vector');
end;
