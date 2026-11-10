function [ys, params, check] = olg80_bc_steadystate(ys, exo, M_, options_)
% Stationary equilibrium of olg80_bc.mod, always with a >= 0.
% The constrained household choices come from lifecycle_bc.m; use
% the steady_state_model block in olg80.mod for the unconstrained life cycle.

params = M_.params;	% Calibrated parameters from the MOD file.
check  = 0;	% Zero indicates a valid stationary solution.
pn     = cellstr(M_.param_names);	% Names used to retrieve parameters.
get    = @(name) params(strcmp(pn, name));	% Look up each calibrated value.
% Parameters read from the calibration
alpha = get('alpha');	% Capital share in production.
beta  = get('beta');	% Annual discount factor.
delta = get('delta');	% Annual capital depreciation.
sigma = get('sigma');	% Relative risk aversion.
Lbar  = get('Lbar');	% Total efficient labor.
J = sum(~cellfun('isempty', regexp(pn, '^z_[0-9]+$', 'once')));	% Number of cohorts.
% Efficiency profile and cohort masses by age
z  = zeros(J,1);	% Efficiency by age.
Nj = zeros(J,1);	% Population share by age.
for j = 1:J
	z(j)  = get(sprintf('z_%d',j));
	Nj(j) = get(sprintf('N_%d',j));
end

% Interest rate that equates asset supply and capital demand
excess = @(r) asset_excess(r, alpha, beta, delta, sigma, Lbar, z, Nj);
[rss, fval, flag] = fzero(excess, [-min(.02,delta/2) .10]);
if flag <= 0 || abs(fval) > 1e-8
	check = 1;
	return
end
[~, kss, yss, wss, css, ass] = asset_excess(rss, alpha, beta, delta, sigma, Lbar, z, Nj);

en = cellstr(M_.endo_names);
idx = @(name) find(strcmp(en, name), 1);
ys = zeros(M_.orig_endo_nbr, 1);
% Fill the steady-state vector in the order of the model variables
ys(idx('Z')) = 1;
ys(idx('y')) = yss;
ys(idx('w')) = wss;
ys(idx('r')) = rss;
ys(idx('k')) = kss;
ys(idx('C')) = Nj'*css;
for j = 1:J
	ys(idx(sprintf('c_%d',j))) = css(j);
end
for j = 1:J-1
	ys(idx(sprintf('a_%d',j))) = ass(j);
end
end

function [gap, k, y, w, c, a] = asset_excess(r, alpha, beta, delta, sigma, Lbar, z, Nj)
% Firms demand capital and pay wages at the candidate interest rate.
k = (alpha*Lbar^(1-alpha)/(r+delta))^(1/(1-alpha));
y = k^alpha*Lbar^(1-alpha);
w = (1-alpha)*y/Lbar;
% Household side: life-cycle problem solved in lifecycle_bc.m
[c, a] = lifecycle_bc(w*z, 1+r, beta, sigma, 0);

% Excess asset supply over capital demand
gap = Nj'*a-k;
end
