function sequences = gibbs_sequences_of_nodes_with_evidence(...
    tree,sequences, F, Q, pi, t, R, X, annealing_prob, read_partition)
% samples the sequences values of cells with reads associated with them
% (and their ancsestors).  Put 0s in all other sequences not tied to
% evidence.

% find subtree with evidence
% T = zeros(size(tree,1),2);
% T(:,1) = tree(:,1);

%  read_partition structure:
%    'raw_from_uniq'
%    'uniq_from_raw'

is_codon = (size(Q,1) > size(R,1));
if is_codon
    global codon2nt;
    codon2ix = [codon2nt(1,:) ; 4+codon2nt(2,:); 8+codon2nt(3,:)]';
    W = repmat(int16(1:64), 3,1)';
    ix = sub2ind([64, 12], W(:), codon2ix(:));
    iy = zeros(64,12);  iy(ix) = 1;
else 
    iy = 1;
end

singles = zeros(size(tree,1), size(Q,1), size(sequences,2)); % size: [M,B,L]
live_cells = find(tree(:,3) == F);
for j=1:length(live_cells)
    uniq_ix = t(read_partition.raw_from_uniq) == live_cells(j);
    rawrep_ix = read_partition.raw_from_uniq(uniq_ix);
    %ix = t == live_cells(j); % check wrt old way
    if ~isempty(uniq_ix)
        y = multi_read_to_counts(X(rawrep_ix,:), 5, read_partition.uniq_weights(uniq_ix));
        %old_y = multi_read_to_counts(X(ix,:), 5); % check wrt old way
        %assert(all(all(y==old_y))); % check wrt old way
        Z = counts_to_potential(y, log(R)); % 4xL
        if is_codon
            Z = reshape(Z,12, []);
        end
        singles(live_cells(j), :, :) = iy*Z;
    end
end
%T = fill_counts(T);
evidence = nodes_with_evidence(tree,t);

% singles = size(sequences,2); % uncomment if you want to sample from prior
% [sequences prob] = tree_sample_for_phylo_unlog(tree(:,1), exp(singles), Q, pi');
% sample the subtree with the evidence
sequences = zeros(size(sequences));
tree_ = tree;
tree_(~evidence,1)=-1;
sequences(evidence,:) = tree_sample_for_phylo_unlog(...
    clean_tree(tree_(:,1)), exp(singles(evidence,:,:)), Q, pi', [], annealing_prob);
% NOTE:  The sequences of the cells that have no evidence are 0!
sequences = int16(sequences);
end

function total = multi_read_to_counts(X, B, weights)    
    
    if ~exist('weights', 'var')
        weights = ones(size(X,1),1);
    end
    
    total = zeros(B, size(X,2));
    for i = 1:size(X,2)
        total(:,i) = accumarray(X(:,i), weights, [B, 1]);
    end
    
    % old way
    if all(weights == 1)
        unweighted_total = histc(X, 1:B, 1);
        assert(all(all(total == unweighted_total)));    	
    end
    
end


function z = counts_to_potential(y, R)
    % R is in log space.
    z = R'*y;  % we assume here that the parents are the *columns* of R
    %z = y*(epsilon.inv-epsilon.log);
end




function sequences = gibbs_sequences_of_dead_cells_old(tree,sequences, F, Q, pi, t, R, X, Qs, rate_class)

low_number = -1000;

global codon2nt;

M = size(tree,1);
B = size(Q,1);
L = size(sequences,2);

% find singleton potentials based on reads.
% singles(k,:,:) is a 64xL potential.
codon2ix = [codon2nt(1,:) ; 4+codon2nt(2,:); 8+codon2nt(3,:)]';
W = repmat(int16(1:64), 3,1)';
ix = sub2ind([64, 12], W(:), codon2ix(:));
iy = zeros(64,12);  iy(ix) = 1;

singles = zeros(M, B, L); 
live_cells = find(tree(:,3) == F);
for j=1:length(live_cells)
    ix = t == live_cells(j);
    if sum(ix) > 0
        y = multi_read_to_counts(X(ix,:), 5);
        Z = counts_to_potential(y, log(R)); % 4xL
        Z = reshape(Z,12, []);
        singles(live_cells(j), :, :) = iy*Z;
    end
end

junk = zeros(B,B*L);

logQ = reshape(log(Q),B,B,L);
logQs = reshape(log(Qs),B,B, []);

% logQ = log(reshape(permute(reshape(Q,B,B,L), [2 1 3]), B, B*L));

% for each node v, in random order
for v=randperm(M)    
    
    Z = zeros(1,B,L);
    
    % find all child sequences
    ix = (tree(:,1) == v);
    if sum(ix) > 0
        % collect transition counts from them
        y = multi_read_to_counts2(sequences(ix,:), B+1);

        % create a BxL singleton potential over the sequences of v        
        for k=1:size(Qs,3)
            iy = rate_class == k;
            Z(1,:,iy) = logQs(:,:,k)'*y(:,iy);
        end
    end

    if tree(v,1) == 0
        sequences(v,:) = tree_sample_for_phylo(0, singles(v,:,:)+Z, junk, pi');
    else
        parent_ix = double(sequences(tree(v,1), :)) + (0:B:(B*L-B));
        sequences(v,:) = tree_sample_for_phylo(0, singles(v,:,:)+Z, junk, Q(:, parent_ix));
    end        
end

%sequences = int16(sequences);


end
   