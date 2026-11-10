function [c,a] = lifecycle_bc(income,R,beta,sigma,initial_assets)
% Deterministic household optimum, constant prices, finite remaining life.
% With a>=0, every prefix budget bounds current consumption. The smallest
% bound determines the first saving spell; repeat from the next age.

income = income(:);		% remaining income stream, one entry per year
N      = numel(income);		% remaining lifetime
c      = zeros(N,1);		% consumption path
a      = zeros(N,1);		% end-of-period assets
g      = (beta*R)^(1/sigma);	% consumption growth factor from the Euler equation
previous = initial_assets;	% assets inherited from the previous age

for j = 1:N
	% Years from the current age
	h = (0:N-j)';

	% Present value of remaining income, discounted at R
	pv = income(j:end)./R.^h;

	% Consumption allowed by each prefix budget (from now to date n)
	bounds = (R*previous+cumsum(pv))./cumsum((g/R).^h);

	% The tightest remaining prefix budget enforces nonnegative assets.
	c(j) = min(bounds);

	% Budget constraint gives end-of-period assets
	a(j)     = R*previous+income(j)-c(j);
	previous = a(j);
end
end
