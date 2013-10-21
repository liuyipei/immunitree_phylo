function [k prob] = many_mnrnd(p, k, should_I_normalize)
% function [k prob] = lmnrnd(l, k)
%
%  Samples from a multi multinomial distribution
%   l is an *already normalized* probability matrix, each *column* a different set of
%     probabilities to sample from
%   k is a vector of the sampled value(s), which can be forced to be what
%     you want.
%   prob = the probability for the outcome in each multinomial

    [K L] = size(p);
    
    % l is a row vector
    if k(1) == 0 
        cdf = cumsum(p,1);
        if nargin < 3 || should_I_normalize == true
            Z = cdf(end,:);
        else
            Z = 1;
            assert( rand>0.01 || abs(cdf(end,1)-1)<1e-5);
        end
        z = rand(1,L).*Z;  
        k = sum(cdf < z(ones(1,K), :), 1) + 1;
        if nargout>1 
            muls = [0:K:(L-1)*K];
            prob = p(muls+k)./Z;
        end
    elseif k(1)>0
        if nargin < 3 || should_I_normalize == true
             Z = sum(p,1);
        else
             Z = 1;             
        end
        muls = [0:K:(L-1)*K];
        prob = p(muls+k)./Z;
    else %if k(1) == -1
        p = p+1e-10*rand(K,L); % resolve ties arbitrarily
        [l_max, k] = max(p,[],1);
        % if you want to normalize
        if nargin < 3 || should_I_normalize == true
            Z = sum(p,1);
        else
            Z = 1;
        end        
        if nargout>1 
            prob = l_max./Z;
        end

    end
end

