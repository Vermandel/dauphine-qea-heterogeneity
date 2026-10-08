function [ys, params, check] = olg80_steadystate(ys, exo, M_, options_)
% Stationary equilibrium of olg80.mod. The scalar root is the interest rate
% that clears the asset market. All parameters are read from the calibration
% written at the top of olg80.mod.

params = M_.params;
check = 0;
pn = cellstr(M_.param_names);
get = @(name) params(strcmp(pn, name));
alpha = get('alpha'); beta = get('beta'); delta = get('delta');
sigma = get('sigma'); Lbar = get('Lbar');
J = sum(~cellfun('isempty', regexp(pn, '^z_[0-9]+$', 'once')));
z = zeros(J,1); Nj = zeros(J,1);
for j = 1:J
    z(j) = get(sprintf('z_%d',j));
    Nj(j) = get(sprintf('N_%d',j));
end

excess = @(r) asset_excess(r, alpha, beta, delta, sigma, Lbar, z, Nj);
[rss, fval, flag] = fzero(excess, [-min(.02,delta/2) .10]);
if flag <= 0 || abs(fval) > 1e-8
    check = 1; return
end
[~, kss, yss, wss, css, ass] = asset_excess(rss, alpha, beta, delta, sigma, Lbar, z, Nj);

en = cellstr(M_.endo_names);
idx = @(name) find(strcmp(en, name), 1);
ys = zeros(M_.orig_endo_nbr, 1);
ys(idx('Z')) = 1;   ys(idx('dep')) = delta;
ys(idx('y')) = yss; ys(idx('w')) = wss; ys(idx('r')) = rss;
ys(idx('k')) = kss; ys(idx('C')) = Nj'*css;
for j = 1:J,   ys(idx(sprintf('c_%d',j))) = css(j); end
for j = 1:J-1, ys(idx(sprintf('a_%d',j))) = ass(j); end
end

function [gap, k, y, w, c, a] = asset_excess(r, alpha, beta, delta, sigma, Lbar, z, Nj)
R = 1+r;
k = (alpha*Lbar^(1-alpha)/(r+delta))^(1/(1-alpha));
y = k^alpha*Lbar^(1-alpha);
w = (1-alpha)*y/Lbar;
g = (beta*R)^(1/sigma);
h = (0:numel(z)-1)';
c = sum(w*z./R.^h)/sum(g.^h./R.^h) * g.^h;   % lifetime budget + Euler growth
a = zeros(numel(z)-1,1);
a(1) = w*z(1)-c(1);
for j = 2:numel(z)-1
    a(j) = R*a(j-1)+w*z(j)-c(j);
end
gap = Nj(1:end-1)'*a-k;
end
