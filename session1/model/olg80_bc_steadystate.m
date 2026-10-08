function [ys, params, check] = olg80_bc_steadystate(ys, exo, M_, options_)
% Stationary equilibrium of olg80_bc.mod, with or without the borrowing limit
% (parameter "constrained"). Household choices come from lifecycle_bc.m.

params = M_.params;
check = 0;
pn = cellstr(M_.param_names);
get = @(name) params(strcmp(pn, name));
alpha = get('alpha'); beta = get('beta'); delta = get('delta');
sigma = get('sigma'); Lbar = get('Lbar'); constrained = get('constrained');
J = sum(~cellfun('isempty', regexp(pn, '^z_[0-9]+$', 'once')));
z = zeros(J,1); Nj = zeros(J,1);
for j = 1:J
    z(j) = get(sprintf('z_%d',j));
    Nj(j) = get(sprintf('N_%d',j));
end

excess = @(r) asset_excess(r, alpha, beta, delta, sigma, Lbar, z, Nj, constrained);
[rss, fval, flag] = fzero(excess, [-min(.02,delta/2) .10]);
if flag <= 0 || abs(fval) > 1e-8
    check = 1; return
end
[~, kss, yss, wss, css, ass] = asset_excess(rss, alpha, beta, delta, sigma, Lbar, z, Nj, constrained);

en = cellstr(M_.endo_names);
idx = @(name) find(strcmp(en, name), 1);
ys = zeros(M_.orig_endo_nbr, 1);
ys(idx('Z')) = 1;   ys(idx('T')) = 0;
ys(idx('y')) = yss; ys(idx('w')) = wss; ys(idx('r')) = rss;
ys(idx('k')) = kss; ys(idx('C')) = Nj'*css;
for j = 1:J,   ys(idx(sprintf('c_%d',j))) = css(j); end
for j = 1:J-1, ys(idx(sprintf('a_%d',j))) = ass(j); end
end

function [gap, k, y, w, c, a] = asset_excess(r, alpha, beta, delta, sigma, Lbar, z, Nj, constrained)
k = (alpha*Lbar^(1-alpha)/(r+delta))^(1/(1-alpha));
y = k^alpha*Lbar^(1-alpha);
w = (1-alpha)*y/Lbar;
[c, a] = lifecycle_bc(w*z, 1+r, beta, sigma, 0, constrained);
gap = Nj'*a-k;
end
