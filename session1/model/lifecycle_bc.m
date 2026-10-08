function [c,a] = lifecycle_bc(income,R,beta,sigma,initial_assets,constrained)
% Deterministic household optimum, constant prices, finite remaining life.
% With a>=0, every prefix budget bounds current consumption. The smallest
% bound determines the first saving spell; repeat from the next age.
income=income(:); N=numel(income); c=zeros(N,1); a=zeros(N,1);
g=(beta*R)^(1/sigma); previous=initial_assets;
for j=1:N
    h=(0:N-j)';
    pv=income(j:end)./R.^h;
    bounds=(R*previous+cumsum(pv))./cumsum((g/R).^h);
    if constrained, c(j)=min(bounds); else, c(j)=bounds(end); end
    a(j)=R*previous+income(j)-c(j);
    previous=a(j);
end
end
