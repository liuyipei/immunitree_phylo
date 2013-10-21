function [k ll best] = get_best_canonic_tree(priors, ll, chain, reads)

if ~exist('reads', 'var'), reads = []; end
priors = init_priors(priors);
st = ceil(size(ll,2)/2); % the first 50% of iterations serve as burn in

% calcualate the chosen sample for each chain.
fprintf('Computing canonical %d trees based on chain:\n', size(ll,2));
for it=st:size(ll,2)
    r = chain(it+1);
    if mod(it,1000) == 0, fprintf('%d ', it); end
    collapse_identical = true;
    a = collapse_edges(convert_phylo_tree_to_mutation_tree(r, collapse_identical, reads));
    T = a.tree2;
    % since tree is much shorter now, reates are resampled and likelihood
    % recomputed for the shorter tree/new rates.
    [~, ll(3,it)] = MH_rates(T, r.F, r.rates, priors, -1);
    [ll(4:5,it),~] = MH_mutation_parameters(T, a.sequences, a.t, r.mut_model, priors);
    
    fprintf('_'); % my progress bar: 50 per line
    if mod(it, 50) == 0 || it == size(ll,2)
        fprintf('\n');
    end
end
ll(:,1:st-1) = -Inf;
[~,k] = max(sum(ll));
best = chain(k+1);
assert(best.iteration == k);
fprintf('Done!\n');
end