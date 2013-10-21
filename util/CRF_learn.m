function [factors reconstruction converged ll] = CRF_learn(X, scopes, nIter, step_size, factors)
global cache
cache = [];
sz = size(X);
if ~exist('step_size', 'var'), step_size = 0.001; end
real_marginals = CRF_marginal_counts(X, scopes);

if ~exist('factors', 'var') || isempty(factors)
    factors = cell(length(scopes),1);
    for i=1:length(factors)-1
        factors{i} = log(real_marginals{i}+1); % can't allow zeros!
    end
    factors{end} = log(real_marginals{end});
end

converged = zeros(length(scopes),1);
ll = zeros(nIter, 1);
% go over every parameter and do one gradient descent step
for it=1:nIter
    if mod(it,10)==0, fprintf('] %d [', it); fprintf('%.2g ',converged); end
    [~, expected_marginals] = CRF_expected_counts(sz, scopes, factors);
    for i=1:length(scopes)-1
        sz_factor = sz(scopes{i});
        ix = myind2sub(sz_factor);
        assert(size(ix,1) == numel(factors{i}));        

        if true %it < nIter - 10
            delta = real_marginals{i} - expected_marginals{i};
            converged(i) = max(delta(:));
            delta = max(min(step_size*delta/(10+it), 1), -1);

%             if max(delta(:)) > converged(i) % we are doing worse than before!
%                 % undo
%                 keyboard;
%                 factors{i} = factors{i} - step_size*old_delta;
%                 %
%             end                
%            if converged(i)<1e-5, continue; end

            factors{i} = factors{i} + delta;
            
%             [~, expected_marginals] = CRF_expected_counts(sz, scopes, factors);
%             new_delta = real_marginals{i} - expected_marginals{i};
%             if max(new_delta(:)) > 2*converged(i) % oh - oh 
%                 iy = find (new_delta(:) > 2*delta(:));
%                 for j=iy % each factor entry 
%                     factors{i}(j) = factors{i}(j) - step_size*delta(j) +  0.1*step_size*delta(j);
%                 end
%             end
%                 
            
        else
            for j=randperm(size(ix,1)) % each factor entry
                [~, expected_marginals] = CRF_expected_counts(sz, scopes, factors);
                delta = real_marginals{i}(j) - expected_marginals{i}(j);
                delta = max(min(step_size*delta/(10+it), 1), -1);
                factors{i}(j) = factors{i}(j) + delta;
            end            
        end
        factors{i} = factors{i}-max(factors{i}(:));
    end
    
    joint_prob = CRF_expected_counts(sz, scopes, factors(1:end-1));
    ll(it) = X(:)' * log(joint_prob(:));
    fprintf('ll = %.2f\n', ll(it));
    
    if all(converged<1e-5), break; end
end

fprintf('\n');
reconstruction = CRF_expected_counts(sz, scopes, factors);
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
[factors reconstruction] = CRF_learn(X , {2, 1, []})


%%
ab = reshape(1:20, 5, 4);
bc = reshape(8:-1:1, 4,2);
X = CRF_expected_counts( [5 4 2], {[1 2], [2 3],[]}, {ab, bc, 1640})
Y = CRF_marginal_counts( X     , {[1 2], [2 3],[]});
[factors reconstruction] = CRF_learn(X, {[1 2], [2 3], []}, 120)
%%
ab = reshape(1:20, 5, 4);
bc = reshape(8:-1:1, 4,2);
X = CRF_expected_counts( [5 4 2], {[1 2], [2 3],[]}, {log(ab), log(bc), log(1640)})
%Y = CRF_marginal_counts( X     , {[1 2], [2 3],[]});
[factors reconstruction] = CRF_learn(X, {[1 2], [2 3], []}, 30140)
end