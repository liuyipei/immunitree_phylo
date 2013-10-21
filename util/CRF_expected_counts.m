% CRF learning
% input:  sz - dimensions of full particle
% factors:  a cell array of factors over subsets of the variables
% scopes:  cell array saying the scope of each factor.  
%
% computes the expected counts for every factor entry under each condition
% scenario.  The last "scopes" entry is the conditioned subset of
% variables.  The last "factors" entry the the total entries drawn from
% each of the conditioned scenario.  The default for that entry 
% is ones(sz(scopes{end})) if not specified.
function [big_hist marginals] = CRF_expected_counts (sz, scopes, factors)

if length(scopes) ~= length(factors)
    assert( length(scopes) - length(factors) == 1 );
%    factors{end+1} = ones([sz(scopes{end}) 1]);
    factors{end+1} = zeros([sz(scopes{end}) 1]); % log space
end

global cache ih
cache_is_on = ~isempty(cache) && isequal(cache{end,1}, scopes) && isequal(cache{end,2}, sz);
if ~cache_is_on, ih = myind2sub(sz); end


big_hist = zeros(sz);
for i=1:length(scopes)-1    
    sz_factor = sz(scopes{i});

    if cache_is_on
        ix = cache{i,1};  
    else
        ix = myind2sub(sz_factor);
    end
    
    assert(size(ix,1) == numel(factors{i}));
    for j=1:size(ix,1) % each factor entry

        if cache_is_on
            iy = cache{i,2}{j};
        else
            iy = true(size(ih,1),1);
            for k=1:size(ix,2)
                iy = iy & (ih(:,scopes{i}(k)) == ix(j, k));
            end
        end
        
%         if i == length(scopes) % last factor!
%             max_hist = max(big_hist(iy));
%             big_hist(iy) = big_hist(iy) - max_hist;
%             big_hist(iy) = big_hist(iy) - log(sum(exp(big_hist(iy))));
% %            big_hist(iy) = big_hist(iy) / sum(big_hist(iy));
%         end
        big_hist(iy) = big_hist(iy) + factors{i}(j);        
    end
end

i=length(scopes);   % last factor
    sz_factor = sz(scopes{i});

    if cache_is_on
        ix = cache{i,1};  
    else
        ix = myind2sub(sz_factor);
    end
    
    assert(size(ix,1) == numel(factors{i}));
    for j=1:size(ix,1) % each factor entry

        if cache_is_on
            iy = cache{i,2}{j};
        else
            iy = true(size(ih,1),1);
            for k=1:size(ix,2)
                iy = iy & (ih(:,scopes{i}(k)) == ix(j, k));
            end
        end
        
        max_hist = max(big_hist(iy));
        big_hist(iy) = big_hist(iy) - max_hist;
        big_hist(iy) = big_hist(iy) - log(sum(exp(big_hist(iy))));
        big_hist(iy) = big_hist(iy) + factors{i}(j);        
    end

big_hist = exp(big_hist);
if nargout > 1
    marginals = CRF_marginal_counts(big_hist, scopes);
end


end

    
    
    



function test()
%%
a = (1:5)';
b = (3:2:9)';
X = CRF_expected_counts( [5 4], {2, 1, []}, {b, a})
%%
X = CRF_expected_counts( [5 4], {2, 1, []}, {b, a, 300})
%%
X = CRF_expected_counts( [5 4], {2, 1, 1}, {b, a})
%%
X = CRF_expected_counts( [5 4], {2, 1, 1}, {b, a, [3 4 3 2 1]'})
%%
X = CRF_expected_counts( [5 4], {2, 1, 2}, {b, a})

X = CRF_expected_counts( [5 4], {2, 1, 2}, {b, a, [4 3 2 1]'})

%%
a_ = [1 2 1 2 1]';
X = CRF_expected_counts( [5 4], {2, 1, 1, []}, {b, a, a_, 100})
X = ceil(X)


%%  3-dimensionsal
ab = reshape(1:20, 5, 4);
bc = reshape(8:-1:1, 4,2);
X = CRF_expected_counts( [5 4 2], {[1 2], [2 3],[]}, {ab, bc, 1640})
Y = CRF_marginal_counts( X     , {[1 2], [2 3],[]});








end