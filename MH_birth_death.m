%function [tree sequences t codes] = MH_birth_death(tree, sequences, F, Q, rates, t, reads, R)
% 
function [tree sequences t codes] = MH_birth_death(tree, sequences, F, Q, rates, t, reads, R)

global greedy % different behavior in "greedy mode"

codes = zeros(1,4);  % some output stats

% double the size of the sequence matrix for faster memory allocation.
% Later on, 'sequences' is allowed to more rows than the number of cells in
% the tree, but only the rows that are up to the num of cells in the tree 
% have any meaning, and at the end of this function we truncate the sequences 
% matrix so it once again matches the size of the tree.    
sequences = [sequences; zeros(size(tree,1), size(sequences,2))];

% number of alive cells (if the parent is "-1" the cell was deleted)
nAlive = sum(tree(:,3) == F & tree(:,1) ~= -1); 

% precompute the set of children for every node
M=size(tree,1);
children = cell(M,1);
for i=1:M
    if tree(i,1) == -1, continue; end
    children{i} = find(tree(:,1) == i);
end

accepted = true;
evidence = nodes_with_evidence(tree,t);

% old version, run on every node
% for v=randperm(size(tree,1)) 
%    if accepted, ix = find(tree(:,1)~= -1); end

% current version:  run K times, where K is the number of nodes tied to 
% evidence at the beginning of this move.
for i=1:sum(evidence)
    % update (if necessary) list of nodes tied to evidence and choose one
    % of them randomly
    if accepted, ix = find(evidence); end     
    v = ix(ceil(rand*length(ix))); 

    % propose a birth-death move centered on node v, compute log likelihood
    % difference and forward transition probability for MH.
    
    % 'code' and 'mv' are input/output parameters that guide the executed 
    % (we feed it to the reverse move so it knows where to come back to)
    code = [v 0 0 0]; % 
    mv = struct('reads_in_u', [], 'child_seq', zeros(1,size(sequences,2)));

    % propose a birth/death move
    [tree_, child_seq, t_, transition_prob, code, mv, ll_diff] = ...
        birth_death(tree, sequences, code, mv, Q, F, rates, t, reads, R, nAlive, children);
    
    % The bias term is due to the fact that 'v' is chosen uniformly over 
    % the number of nodes, and that number changes slightly with a death/birth.
    nAlive_ = nAlive;
    children_old = children{v};
    if ~isempty(child_seq) % a new cell was born, matching the last row in 'tree_'        
        child_id = size(tree_,1);
        sequences(child_id,:) = child_seq;
        bias = length(ix)/(length(ix)+code(4)); % code(4) says whether new cell has reads
        if tree_(end,3) == F, nAlive_ = nAlive+1; end
        children{v} = [children{v}; child_id];
        children{child_id} = [];        
    else % a cell was deleted
        bias = length(ix)/(length(ix)-evidence(code(4)));
        if code(3) == F, nAlive_ = nAlive-1; end
        children{v} = find(tree_(:,1) == v);
    end
    
    % compute the reverse transition probability and reverse likelihood
    % difference. 
    % note: ll_diff is only computed upon a 'birth' move, and it is set to 
    % 0 by a 'death' move, and hence ll_diff-rev_ll_diff is the actual 
    % likelihood ratio (in log scale)
    [~,~,~, rev_transition_prob, ~, ~, rev_ll_diff] = ...
        birth_death(tree_, sequences, code, mv, Q, F, rates, t_, reads, R, nAlive_, children);

    % compute MH acceptance probability, including the bias term
    acceptance_prob = bias*exp(ll_diff-rev_ll_diff)*rev_transition_prob/transition_prob;
       
    accepted = rand<acceptance_prob || (greedy == -1 && acceptance_prob > 0.5);
    if accepted % accept changes to tree
        tree = tree_;
        t = t_;
        nAlive = nAlive_;

        if code(2) == -1  % a birth
            evidence(child_id) = code(4); % code(4) says whether new cell has reads
            codes(1) = codes(1)+1;
            if code(3) == 1, codes(4) = codes(4)+1; end

        else  % a death
            evidence(code(4)) = false;
            codes(2) = codes(2)+1;
        end
        
    else % rejected
        children{v} = children_old;
        codes(3) = codes(3)+1; % rejection count
    end
end

% truncate 'sequences' matrix
sequences = sequences(1:size(tree,1), :); 

% remove deleted nodes from tree
[tree sequences t] = clean_tree(tree, sequences, t); 

end




% this is the actual move. 
% 'code' and 'mv' as inputs are the outputs of the "forward move"
% they are empty (except code(1)) if we are just constructing the forward move.  
% code(1) = the active node, 'v'
% if move was birth: set code(2) = -1 
%                        code(3) = (is new child dead?)
%                        code(4) = (does new child have reads?)
% if move was death: set code(2:3) = birth and death time of deleted node
%                        code(4) = id of deleted node
% if code(2) = 0, it means there was no "forward move" yet, we 
% are making it right now.
%
% v is the active cell.  
% The move deletes one of its *leaf* children or create a new (leaf) child
% if v is alive:
%    in case of creation: some of v's reads might move to the new child.
%    in case of deletion: all of the deleted leaf's reads will move back to v.
% if v is dead:
%   we are not allowed to delete a leaf that has reads (because reads
%   cannot be assigned to dead cells)
%
% mv has two fields, filled only when deleting a leaf
%     mv.reads_in_u = ids of reads in the deleted leaf
%     mv.child_seq = sequence of deleted leaf
function [tree child_seq t transition_prob code mv ll_diff] = ...
    birth_death(tree, sequences, code, mv, Q, F, rates, t, reads, R, nAlive, children)

v = code(1);
child_seq = [];

transition_prob = 1;

% Create a list of leaf children that are allowed to be deleted
ix = children{v};  % all children of v
cand = false(1, length(ix));
for i=1:length(ix)    
    % If v is alive include all leaf children.
    % If v is dead, include all leaf children with no reads.
    u = ix(i); 
    if isempty(children{u}) && (tree(v,3) == F || isempty(find(t == u, 1)))
        cand(i) = 1;
    end
end
ix = ix(cand); % ix now has all valid candidates for deletion


% if forward move was "birth" then confirm reverse move has children.
assert(~isempty(ix) || code(2) ~= -1); 

if ~isempty(ix) % at least one candidate for deleteion
    transition_prob = 0.5; % chance for deletion
    
    % if this is the forward move, decide (by sampling) if we go for deletion.
    % if this is the reverse move, go for deletion only if forward move was creation.
    if (code(2) == 0 && rand < transition_prob) || code(2) == -1
        % erase empty node 
        u = ix(ceil(length(ix)*rand));  % uniformly choose a node to delete 
        transition_prob = transition_prob*1/length(ix);
       
        tree(u,1) = -1; % delete u
        code(2:3) = tree(u,2:3); % record birth/death times for reverse
        code(4) = u; 

        was_dead = (tree(u,3) < F); % was deleted node dead?

        ll_diff = 0;  % not computed for deletion moves
        
        % it's more likely to have less live cells
        if ~was_dead, 
            % record the reads that belong to u and u's sequence
            mv.reads_in_u = find(t==u);
            mv.child_seq = sequences(u,:);
            
            % Associate reads of u to the parent v.        
            t(mv.reads_in_u) = v;
        end
        return;
    end
    % we are going for creation!  
    transition_prob = 1-transition_prob;    
end

% Prepare a proposal to add an empty child:

ll_diff = 0;

% sample birth time
birth_time = code(2);
constraint = tree(v,3)-tree(v,2);
if birth_time == 0
    birth_time = tree(v,2) + mod(exprnd(1/rates.birth), constraint); % conditional probability
end

% sample death time (given birth time)
death_time = code(3);
if death_time == 0
    death_time = birth_time + exprnd(1/rates.death);
end
death_time = min(death_time,F); 

% sample sequence of child from prior
B = size(Q,1); L = size(sequences,2);
assert(size(Q,2) == B*L);
parent_ix = double(sequences(v,:)) + (0:B:(B*L-B));
parent_seq = sequences(v,:);
child_seq = many_mnrnd(Q(:,parent_ix),double(mv.child_seq),false);
if all(mv.child_seq == 0)
    while all(child_seq==parent_seq)
        child_seq = many_mnrnd(Q(:,parent_ix),double(mv.child_seq),false);
    end
end
child_seq = many_mnrnd(Q(:,parent_ix),double(mv.child_seq),false);

child_seq = int16(child_seq);

% note: because we sampled the child sequence from the prior, no need to 
%       add this sampling probability to either the transition probability 
%       or the likelihood difference, because both additions cancel each 
%       other out.
%
%       Other terms cancel out as well, due to the similarity between the
%       prior and the way we sample our birth move

%TODO: allocate in advance space for the tree - it takes a lot of time to append.
tree(end+1,:) = [v birth_time death_time];

ll_diff = ll_diff - (rates.birth+rates.death)*(death_time-birth_time);

transition_prob = transition_prob ...
    *exp(-rates.birth*(birth_time-tree(v,2)))/(1-exp(-rates.birth*constraint)) ....
    *exp(-rates.death*(death_time-birth_time));

code(2) = -1;  % mark a "birth" event

% if we created a live cell, move some of v's reads to it
 if death_time ==  F,     
    % read assignment likelihood: it's more likely to have less live cells
    ll_diff = ll_diff + length(t)*log(nAlive/(nAlive+1));
    code(3) = 1;
    
    % if parent has reads, move some of them to child
    reads_in_v = find(t==v);
    if isempty(reads_in_v), assert(isempty(mv.reads_in_u)); return; end
    
    seqs = [sequences(v,:); child_seq];
    if size(Q,1) > size(R,1) % AA model
        seqs = codons2seqs(seqs);
    end
    
    % compute the log likelihood weights, for every read, between being
    % associated with v or the new child 
    mut = (seqs(1,:) ~= seqs(2,:));
    logR = log(R); % R is the mutation matrix
    LLs = zeros(2, length(reads_in_v));
    if sum(mut) > 0
        for i=1:length(reads_in_v)
            for j=1:2
                LLs(j,i) = sum(logR(5*(seqs(j,mut)-1)+ reads(reads_in_v(i),mut)));
            end
        end
    end
    
    % if this is a reverse move, we know which reads must move to the child
    % (they are listed in mv.reads_in_u).  Mark in the array 'force' '1' for
    % every read of v that stays in v, and '2' for the reads that move.
    forced = zeros(1,length(t));
    assert(code(3) ~= 0) % Yi thinks this must be the case -- Aug 2012
    if code(3) ~= 0 % this is a reverse move    
        forced(reads_in_v) = 1;
        forced(mv.reads_in_u) = 2;    
    end
    forced = forced(reads_in_v);

    % If this is the forward move, 'force' will have zeros, which tells the
    % 'many_lmnrnd' to sample t_ according to the weights in LLs.    
    [t_ prob] = many_lmnrnd(LLs, forced);
    t(reads_in_v(t_==2)) = size(tree,1); % move the reads where t_=2 to child
    code(4) = ~isempty(find(t_==2, 1)); % at least one read moved
        
    % change ll_diff and transition_prob accordingly.
    ll_diff = ll_diff + sum([-1 1]*LLs(:,t_==2));
    transition_prob = transition_prob * prod(prob);    
    
 else 
    code(3) = 0;
 end

end
