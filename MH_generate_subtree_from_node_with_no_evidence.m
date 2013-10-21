function [tree sequences t codes] = MH_generate_subtree_from_node_with_no_evidence(tree, ...
    sequences, t, rates, F, Q, pi)


B = size(Q,1); 
L = size(sequences,2); 
N = length(t);
assert(size(Q,2) == B*L);

evidence = nodes_with_evidence(tree,t);
% roots are non-evidence nodes whose parent are evidence nodes.
if length(evidence)>1
    roots = [false; (~evidence(2:end) & evidence(tree(2:end,1))) ];
else
    roots = false; 
end

roots = find(roots);

alive = tree(:,3)==F;
T = [tree(:,1) alive];
T = fill_counts(T);
M = sum(alive);
M0 = M;
C0 = size(tree,1);

S = [tree(:,1) zeros(size(tree,1),1)];
S(roots,2) = roots;
S = fill_down(S);

sequences = [sequences; zeros(size(sequences), 'int16')];
nSeqBuff = size(sequences,1);

codes = zeros(1,2);
for r=randperm(length(roots))
    % count live descendents of that root.
    nCells = sum(S(:,3)==roots(r));
    nAlive = T(roots(r),3);    
    
    F_ = F-tree(roots(r),2);
    %fprintf('%d %d %.4f\n', r, roots(r), F_);
    parent_ix = double(sequences(tree(roots(r),1),:)) + (0:B:(B*L-B));
    maxNodes = nCells+10;
    [tree_ sequences_ F__]  = generate_sequences_from_prior(maxNodes, F_, rates, Q, Q(:,parent_ix));
    if (F__<F_), codes(2) = codes(2)+1; continue; end % too many new live cells
    
    nCells_ = size(tree_,1);
    nAlive_ = sum(tree_(:,3)==F_);    
%    fprintf('%d more alive, %d more total.  ', nAlive_-nAlive, nCells_-nCells);
    if log(rand) < N*(log(M)-log(M-nAlive+nAlive_))       
        tree_(:,2:3) = tree_(:,2:3)+tree(roots(r),2);
        tree_(:,1) = tree_(:,1) + size(tree,1);
        tree_(1,1) = tree(roots(r),1);

        % erase old nodes
        tree(S(:,3)==roots(r),:) = -1;        

        if size(tree_,1) == 1 % no new nodes
            tree(roots(r),:) = tree_;
            sequences(roots(r),:) = sequences_;
        else         % append new nodes
            while nSeqBuff-size(tree,1) < size(tree_,1)
                sequences = [sequences; zeros(size(sequences), 'int16')];
                nSeqBuff = size(sequences,1);
            end
            sequences((size(tree,1)+1):(size(tree,1)+size(tree_,1)),:) = sequences_;
            tree = [tree; tree_];
        end
        
        % sanity check: number of alive nodes is consistent
        %assert(sum(tree(:,3)==F) == M-nAlive+nAlive_);
        M = M-nAlive+nAlive_;
%        fprintf('accepted\n');
        codes(1) = codes(1)+1;
    else 
%        fprintf('rejected\n');
        codes(2) = codes(2)+1;
    end
end
sequences = sequences(1:size(tree,1),:);
[tree sequences t] = clean_tree(tree, sequences, t);

% sample (from the prior) the sequences for the rejected subtrees
% Those sequences are those with 0 value.
sequences = tree_sample_for_phylo_unlog(tree(:,1), size(sequences,2), Q, pi', sequences);

fprintf('%.2f%% accepted, %.2f%% rejected, %d more cells, %d more alive\n', ...
    100*codes(1)/length(roots), 100*codes(2)/length(roots), size(tree,1)-C0, M-M0);

end



