function [tree total_ll_diff codes] = MH_auxiliary_graph(tree, sequences, Q, t) 
% So what happens here?
% 1) calculate hamming distance between all cells to all cells
% 2) for each cell v, construct a pool of candidate parents. a cell is a 
%    candiate parent with prob proportional to exp(-k*distance), and if it 
%    alive when v is born.  Force v's current parent to be included.
% 3) Choose v's parent uniformally from the pool.  For faster mixing, do
%    not consider v's current parent.
% 4) denote by ix the cells in the pool.  denote by x the current state.
%    ix is an auxiliary variable.  It does not change in the move.
%    the acceptance prob: P(x')/P(x) * P(ix|x')/P(ix|x) * T(x'->x)/T(x->x')
%    P(x')/P(x) is calculated based on the different mutations creating v.
%    T(x'->x)/T(x->x') is 1 (both are 1/[length(ix)-1] )
%    P(ix|x')/P(ix|x) = P(v's current parent is in the pool)/P(v's new parent is in the pool)
%         [this is because we forced the current parent to be in ix]

global greedy

assert(size(tree,1) == size(sequences,1));
total_ll_diff = 0;
codes = zeros(1,3);
M = size(tree,1);
B = size(Q,1); 


% calculate hamming distance between all nodes to all nodes
[seqs I J] = unique(sequences, 'rows');
h = hist(J, 1:length(I));
%dist = pdist(double(seqs), 'hamming')*size(seqs,2); % pdist is slow

dist = matmult_as_pairwise_hamming(seqs);
if isempty(dist), fprintf('all the sequences are the same!\n'); return; end

% generate an undirected graph of nodes, chance of a node is proportional 
% to the distance.  

% the expected number of edges is:
mean_degree = 35;
expected_num_of_edges = mean_degree*M/2;

dist = exp(-2*dist);
%probs = squareform(dist)+eye(length(I)); % used to use pdist/hamming,
%which was not in square form
%assert(max(max(squareform(dist)-dist2))<100*eps); % old correctness check
probs = dist+eye(length(I));
probs = 0.5*expected_num_of_edges * probs/(h*probs*h'-length(J));


% makes sure we won't move nodes as children of deleted nodes
erased = (tree(:,1) == -1);
assert(sum(erased)==0);
% probs(erased,:) = 0;
% probs(:,erased) = 0;

%%% Since lifespans don't change, we compute them in advance
DG = false(M,M);
times = [tree(:,2) (1:M)'; tree(:,3) -(1:M)'];
times = sortrows(times);
buff = false(M,1);
for i=1:length(times)
    if times(i,2)<0
        buff(-times(i,2)) = false;
    else
        DG(:,times(i,2)) = buff;
        buff(times(i,2)) = true;
    end
end
%%% 

%%% Only consider detaching nodes that are tied to evidence
[evidence evidence_quantity]= nodes_with_evidence(tree,t); 
direct_evidence = false(size(tree,1),1); % precompute direct evidence for performance
direct_evidence(t) = true;

accepted = true;
for i=1:sum(evidence) %
    if accepted, pool = find(evidence); end
    v = pool(ceil(rand*length(pool))); 
%%% 
%for v=randperm(size(tree,1))

    if tree(v,1) <= 0, continue; end % ignore root and deleted nodes
    code = [v 0];

    % for each node, choose a neighbour uniformly at random, and move the node
    % there
    [tree_ code] = auxiliary_graph(tree, code, probs, J, DG);
    if code(2) == 0, codes(3) = codes(3)+1; continue; end
    
    ll_diff = 0;
    a = tree(v,1);
    b = tree_(v,1);
    assert(b~=a);
    ix = find(sequences(b,:) ~= sequences(a,:));
    if ~isempty(ix)
        ll_diff = ll_diff + ...
            -sum(log(Q(mysub2ind(B, sequences(v,ix), sequences(a,ix), ix)))) ...
            +sum(log(Q(mysub2ind(B, sequences(v,ix), sequences(b,ix), ix))));
    end

    bias = 1;
    
%%% update evidence structure, and calculate bias

%%% perf optimization/simplification: new strat: get all ancestors of old parent(a) and
%%% new parent(b). Then, update the net evidence with each.
    p = b; ev_ft = [];
    while p > 0
        ev_ft = [ev_ft p]; 
        assert(p ~= tree(p,1)); % no fixed points in traversing up the tree!
        p = tree_(p,1);
    end    
    
    p = a; ev_tf = [];    
    while p > 0
        ev_tf = [ev_tf p]; 
        assert(p ~= tree(p,1));
        p = tree_(p,1); 
    end
    
    evidence_quantity(ev_ft) = evidence_quantity(ev_ft) + evidence_quantity(v); 
    evidence_quantity(ev_tf) = evidence_quantity(ev_tf) - evidence_quantity(v);     
    evidence(ev_ft) = evidence_quantity(ev_ft) > 0;
    evidence(ev_tf) = evidence_quantity(ev_tf) > 0;
    bias = length(pool)/(length(pool)-length(ev_tf)+length(ev_ft));
%%% 
    
    acceptance_prob = bias*exp(ll_diff)*probs(J(v),J(a))/probs(J(v),J(b));
%    fprintf('nSeqDiff: %d  ll_diff = %.2f trans_diff = %.2f  acc=%.2f\n', length(ix), ll_diff, log(probs(J(v),J(a))/probs(J(v),J(b))), acceptance_prob);
    accepted = rand<acceptance_prob;
    if accepted || (greedy == -1 && ll_diff > 0)
        codes(1) = codes(1)+1;
        tree = tree_;
        total_ll_diff = total_ll_diff + ll_diff;
    else
        codes(2) = codes(2)+1;
        %%% revert changes in evidence
        evidence_quantity(ev_ft) = evidence_quantity(ev_ft) - evidence_quantity(v); 
        evidence_quantity(ev_tf) = evidence_quantity(ev_tf) + evidence_quantity(v);     
        evidence(ev_ft) = evidence_quantity(ev_ft) > 0;
        evidence(ev_tf) = evidence_quantity(ev_tf) > 0;
        %%%
    end
end
assert(isequal(evidence,nodes_with_evidence(tree,t))); % is evidence still consistent

end

function [tree code] = auxiliary_graph(tree, code, probs, J, DG)

M = size(tree,1);
v = code(1);
assert(code(2) == 0);

% DONE: Switch order:  first filter the columns of G, and then randomize.
% choose a pool of parents whose lifespan overlap v's birth
ix = DG(:,v);
% ix_ = ( (tree(:,2)<tree(v,2)) & (tree(v,2)<tree(:,3)) );
% assert(isequal(ix,ix_));
ix(tree(v,1)) = false;
ix = find(ix);
G = (rand(1,length(ix))<probs(J(v),J(ix)));
ix = ix(G);

u = ceil(length(ix)*rand);
if u == 0,  return; end
code(2) = tree(v,1);
tree(v,1) = ix(u);
%sanity(tree);
end


function sanity(tree)
    for u=2:size(tree,1);
        assert(tree(u,2)>tree(tree(u,1),2));
        assert(tree(u,2)<=tree(tree(u,1),3));
        assert(tree(u,2)<=tree(u,3));
    end
end

function [pairwise_hamming] = ...
    matmult_as_pairwise_hamming(rows2compare) % the alphabet consists of 1..B

    B = int16(max(max(rows2compare)));
    augmented_matrix = zeros(size(rows2compare,1),B*size(rows2compare,2));
    muls = int16(0:(size(rows2compare,2)-1))*int16(B);
    for i = 1:size(rows2compare,1);
        augmented_matrix(i,muls+int16(rows2compare(i,:))) = 1;
    end
    pairwise_hamming = size(rows2compare, 2) - (augmented_matrix * (augmented_matrix'));
end
