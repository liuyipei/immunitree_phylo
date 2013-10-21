% input:  X = a histogram of counts over 2 variabels
% output: p = p-value for X drawn from a distribution in which the 2
%             variables are indpendent. small p <--> not independent
function [p stats] = chi2independence(X)

% x <- distribution over non-zero columns of X
epsilon = 0;%0.00001;
X = X+epsilon;
N = sum(X(:));
x = sum(X,1)/N;
X(:,x==0) = [];
x(x==0) = [];

% y <- distribution over non-zero rows of X
y = sum(X,2)/N;
X(y==0,:) = [];
y(y==0) = [];

% expected counts according to independence assumption
expected = y*x*N; 
X = X-epsilon;

% goodness of fit comparing actual counts (X(:)) with expected(:)
bins = (1:numel(X))';
[h,p,stats] = chi2gof(bins, 'frequency', X(:),  'ctrs', bins, 'expected', ...
    expected(:), 'nparam', length(x) + length(y) - 2, 'emin', 0);

% alaternative way to construct p value
% return 
p_ = p;
d = [];
for i=1:size(X,1)
    for j=1:size(X,2)
        d = [d; repmat([i j], X(i,j),1)];
    end
end
[~,stats,p] = crosstab(d(:,1), d(:,2));
assert(abs(p-p_)<1e-5);
end


%%
function test()
%%
    % draw roughly 400 counts from a random distribution over 2 variables
    Y = ceil(20*rand(4,10));
    
    % should not be significant
    p_significant = chi2independence(Y)
    
    a = ceil(100*rand(4,1));
    b = ceil(100*rand(10,1));
    mult = a*b'/sum(a)/sum(b);
    X = mult*1000;
    Z = ceil(X)+Y;
    p_not_significant = chi2independence(Z)

    
end