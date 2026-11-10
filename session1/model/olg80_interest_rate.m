function rss = olg80_interest_rate(alpha,beta,delta,sigma,Lbar)
% Solve only the scalar capital-market equilibrium. The life-cycle c_1
% at each trial rate is the same analytic budget/Euler formula shown in
% olg80.mod; the other steady-state variables are built there, not here.
global M_
pn = cellstr(M_.param_names);
J = sum(~cellfun('isempty',regexp(pn,'^z_[0-9]+$','once')));
get = @(name) M_.params(strcmp(pn,name));
z = arrayfun(@(j) get(sprintf('z_%d',j)),(1:J)');
Nj = arrayfun(@(j) get(sprintf('N_%d',j)),(1:J)');
excess = @(r) asset_excess(r,alpha,beta,delta,sigma,Lbar,z,Nj);
[rss,fval,flag] = fzero(excess,[-min(.02,delta/2) .10]);
assert(flag > 0 && abs(fval) < 1e-8, 'Capital market did not clear');
end

function gap = asset_excess(r,alpha,beta,delta,sigma,Lbar,z,Nj)
% Trial factor prices and the analytic life-cycle budget.
R = 1+r;
k = (alpha*Lbar^(1-alpha)/(r+delta))^(1/(1-alpha));
w = (1-alpha)*k^alpha*Lbar^(-alpha);
g = (beta*R)^(1/sigma);
h = (0:numel(z)-1)';
c1 = w*sum(z./R.^h)/sum((g/R).^h);
c = c1*g.^h;
a = zeros(numel(z)-1,1);
a(1) = w*z(1)-c(1);
for j = 2:numel(a)
    a(j) = R*a(j-1)+w*z(j)-c(j);
end
gap = Nj(1:end-1)'*a-k;
end
