function [k prob] = many_lmnrnd(logp, k)
% function [k prob] = lmnrnd(l, k)
%
%  Samples from a multi multinomial distribution
%   l is a log probability matrix, each *column* a different set of
%     probabilities to sample from
%   k is a vector of the sampled value(s), which can be forced to be what
%     you want.
%   prob = the probability for the outcome in each multinomial

    [K L] = size(logp);
    
    % l is a row vector  (this case just takes the max? comment Yi)
    if k(1) == -1
        logp = logp+1e-10*rand(size(logp)); % resolve ties arbitrarily
        [l_max m] = max(logp,[],1);
        k = m;
    else
        l_max = max(logp,[],1);        
    end    
  
    
    p = exp(logp-l_max(ones(1,K), :));

    % normalize to sum to 1 (and correct numerical issues)
    cdf = cumsum(p,1);

    if k(1) == 0  % samples an integer, per column
        % cumsum and find min
        z = rand(1,L).*cdf(end, :);
        k = sum(cdf < z(ones(1,K), :), 1) + 1;
    end
    if nargout>1 
        muls = [0:K:(L-1)*K];
        prob = p(muls+double(k))./cdf(end, :); 
    end
end

