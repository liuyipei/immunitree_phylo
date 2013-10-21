function [k prob] = lmnrnd(logp, k)
% function k = lmnrnd(l, num)
%
%  Samples from a multinomial distribution
%   l is a log probability row vector
%    if nargin<2, k = 0; end

    [l_max, m] = max(logp);
    
    p = exp(logp-l_max);
    % normalize to sum to 1 (and correct numerical issues)
    
    cdf = cumsum(p);       
    
%    k = max([(k==0)*sum(cdf < rand*cdf(end))+1 -m*k k]);
     % cumsum and find min
     if k == 0
         k = sum(cdf < rand*cdf(end)) + 1; 
     elseif k == -1 
         k = m; 
     end
    prob = p(k)/cdf(end);
end


% logp = [log_a log_b]
% assume l_max = max(logp) = log_a
% ==> p = exp( [0 log(b/a) ]) = [1 b/a]
% ==> sum_p = 1+b/a
% ==> norm_p = [a/(a+b) b/(a+b)] 