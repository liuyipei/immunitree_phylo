function [k prob] = super_lmnrnd(logp, num, k)
% function k = super_lmnrnd(l, num, k)
%
%  Samples from a multinomial distribution
%   l is a log probability vector
%   num is the number of samples to draw (default = 1)
%   k force function to choose entry k (over-rule random)

    if nargin < 2, num = 1; end
    if num == 0, k = []; return; end 

    % l is a row vector    
    [l_max, m] = max(logp);
%    k = m; return; % in case you want to be greedy
    
    p = exp(logp-l_max);
    % normalize to sum to 1 (and correct numerical issues)
    sum_p = sum(p);
    norm_p = p/sum_p;
    if nargin == 3 && k > 0, prob = norm_p(k); return; end

    
    % cumsum and find min
    if num == 1
        D = rand < cumsum(norm_p);    
    else
        Z = rand(num,1);
        W = cumsum(norm_p);
        D = (Z(:,ones(1,length(norm_p))) < W(ones(num,1),:) )';
%        D = (repmat(rand(num,1),1,length(norm_p)) < repmat(cumsum(norm_p),num,1) )';    
    end
    k = (length(norm_p)+1) - sum(D);
    prob = norm_p(k);
end
