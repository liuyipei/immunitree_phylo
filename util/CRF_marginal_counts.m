function factors = CRF_marginal_counts(big_hist, scopes)

    sz = size(big_hist);
    global cache ih  
    cache_is_on = ~isempty(cache) && isequal(cache{end,1}, scopes) && isequal(cache{end,2}, sz);
    if ~cache_is_on
        cache = cell(length(scopes)+1,2);
        cache{end,1} = scopes;
        cache{end,2} = sz;
        ih = myind2sub(sz);
    end
    
    
    factors = cell(length(scopes),1);
    for i=1:length(scopes)
        sz_factor = sz(scopes{i});        
        factors{i} = zeros([sz_factor 1]);        

        if cache_is_on
            ix = cache{i,1};
        else
            ix = myind2sub(sz_factor);
            cache{i,1} = ix;
            cache{i,2} = cell(size(ix,1),1);
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
                cache{i,2}{j} = iy;
            end
            
            factors{i}(j) = sum(big_hist(iy));
        end
        %factors{i} = log(factors{i});
    end
end


function test()
%%
a = (1:5)';
b = (2:2:8)';
X = [2     4     6     8;
     4     8    12    16;
     6    12    18    24;
     8    16    24    32;
    10    20    30    40];
factors = CRF_marginal_counts(X , {2, 1, 1})
end