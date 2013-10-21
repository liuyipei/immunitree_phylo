function b = convert_phylo_tree_to_mutation_tree(a, collapse_identical, reads)

if ~exist('collapse_identical', 'var'), collapse_identical = false; end
if ~exist('reads', 'var'), reads = []; end


if collapse_identical
    assert(all((1:size(a.tree,1))' ~= a.tree(:,1)))
    max_collapse_iterations = 100;
    for i = 1:max_collapse_iterations
        old_tree_size = sum(a.t(:,1) > 0);
        [a.tree, a.sequences, a.t, ~, ~] = ... % added 08/2012
            tree_collapse_identical_nodes(a.tree, a.sequences, a.t, [], a.F);
        if sum(a.t(:,1) > 0) == old_tree_size
            % nothing has changed -- hence converged
            break
        end

        if ~isempty(reads)
            Qs = get_Qs(a.mut_model);
            annealing_prob = 1;            
            a.sequences = gibbs_sequences_of_nodes_with_evidence(a.tree, ...
                a.sequences, a.F, Qs, a.mut_model.pi, a.t, a.R, reads, annealing_prob); %            
        else
            break % no convergence issue -- break
        end
    end
    if i == max_collapse_iterations
        fprintf('convert_phylo_tree_to_mutation_tree: max_collapse_iterations (%d) reached!', max_collapse_iterations);
    end
end
assert(all((1:size(a.tree,1))' ~= a.tree(:,1)))
death = a.tree(:,3);

parent = a.tree(:,1);
sequences = a.sequences;
count = hist(a.t,1:size(a.tree,1))'; 
t = a.t;

for k=1:length(parent)
    if parent(k) ~= 0
        if all(sequences(k,:) == sequences(parent(k),:))
            % erase k
            count(parent(k)) = count(parent(k)) + count(k);
            death(parent(k)) = max(death(parent(k)), death(k)); % added
            if ~isempty(t), t(t==k) = parent(k); end
            parent(parent == k) = parent(k);
            parent(k) = -1;  
        end
    end
end
parent2 = [parent(:,1) a.tree(:,2) death];
parent2 = clean_tree(parent2);

parent = [parent count];
[parent sequences t] = clean_tree(parent, sequences, t);
parent = fill_counts(parent);

% erase sub-trees with no observations
ix = find(parent(:,3)==0);
%fprintf('%d nodes in the tree belong to subtrees with no observations.\n', length(ix));
parent(ix,1) = -1;  

% added July 12th 2011.
parent2 = [parent(:,1) parent2(:,2:3)];
parent2 = clean_tree(parent2);
b.tree2 = parent2;

[parent sequences t] = clean_tree(parent,sequences, t);

b.tree = parent;
b.sequences = sequences;
b.t = t;
end


