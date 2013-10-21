function [tree sequences t codes] = MH_birth_death(tree, sequences, F, Q, birth_rate, death_rate, t, reads, R)

global greedy

%total_ll_diff = 0;
codes = zeros(1,4);

sequences = [sequences; zeros(size(tree,1), size(sequences,2))];
nAlive = sum(tree(:,3) == F & tree(:,1) ~= -1);

%%% TODO: confirm that using the sparse thing isn't faster.
M=size(tree,1);
% DG = sparse(1+tree(:,1), 2:M+1, true, M+1, M+1);
% DG = DG(2:end, 2:end);
children = cell(M,1);
for i=1:M
    if tree(i,1) == -1, continue; end
%    children{i} = DG(i,:); 
    children{i} = find(tree(:,1) == i);
end
accepted = true;
evidence = nodes_with_evidence(tree,t);
%%%

%for v=randperm(size(tree,1))
% only consider nodes that are tied with evidence
for i=1:sum(evidence)
%    if accepted, ix = find(tree(:,1)~= -1); end
    if accepted, ix = find(evidence); end

    v = ix(ceil(rand*length(ix)));
    code = [v 0 0 0]; 
    mv = struct('reads_in_u', [], 'child_seq', zeros(1,size(sequences,2)));
    [tree_ child_seq t_ transition_prob code mv ll_diff] = ...
        birth_death(tree, sequences, code, mv, Q, F, birth_rate, death_rate, t, reads, R, nAlive, children);

    if code(2) == 0, assert(false); codes(3) = codes(3) + 1; continue; end       
    % sequences is allowed to grow beyond the size of the tree, but only
    % the sequences that are up to the size of the tree have any meaning,
    % and at the end of this function we trim the sequences array to
    % reflect those of the tree.
    
    % The bias term is due to the fact that choosing 'v' is chosen
    % uniformly over the the number of nodes, and that number changes
    % slightly with a death/birth.
    nAlive_ = nAlive;
    children_old = children{v};
    if ~isempty(child_seq) % a new cell was born, matching the last row in 'tree'        
        child_id = size(tree_,1);
        sequences(child_id,:) = child_seq;
%        bias = length(ix)/(length(ix)+1);
        bias = length(ix)/(length(ix)+code(4)); % code(4) says whether new cell has reads
        if tree_(end,3) == F, nAlive_ = nAlive+1; end
        children{v} = [children{v}; child_id];
        children{child_id} = [];        
    else % a cell was deleted
%        bias = length(ix)/(length(ix)-1);
        bias = length(ix)/(length(ix)-evidence(code(4)));
        if code(3) == F, nAlive_ = nAlive-1; end
        children{v} = find(tree_(:,1) == v);
    end
    
    [~,~,~, rev_transition_prob, ~, ~, rev_ll_diff] = ...
        birth_death(tree_, sequences, code, mv, Q, F, birth_rate, death_rate, t_, reads, R, nAlive_, children);
    acceptance_prob = bias*exp(ll_diff-rev_ll_diff)*rev_transition_prob/transition_prob;
       
    accepted = rand<acceptance_prob || (greedy == -1 && acceptance_prob > 0.5);
    if accepted
        tree = tree_;
        t = t_;
        nAlive = nAlive_;
        if code(2) == -1  % a birth
            evidence(child_id) = code(4); % code(4) says whether new cell has reads
            codes(1) = codes(1)+1;
            if code(3) == 1, codes(4) = codes(4)+1; end
        else              % a death
            evidence(code(4)) = false;
            codes(2) = codes(2)+1;
        end
    else
        children{v} = children_old;
        codes(3) = codes(3)+1;
    end
end

sequences = sequences(1:size(tree,1), :);
assert(size(tree,1) == size(sequences,1));
assert(isequal(evidence,nodes_with_evidence(tree,t)));
[tree sequences t] = clean_tree(tree, sequences, t);

end




function [tree child_seq t transition_prob code mv ll_diff] = ...
    birth_death(tree, sequences, code, mv, Q, F, birth_rate, death_rate, t, reads, R, nAlive, children)

v = code(1);
child_seq = [];

transition_prob = 1;

% nAlive_ = sum(tree(:,3) == F & tree(:,1) ~= -1);
% assert(nAlive == nAlive_);

% collect children that are leaves
% prepare proposals where each one of those children does not exist
ix = children{v};  % all children of v
% ix_ = find(tree(:,1) == v);  % all children of v
% assert(isequal(ix(:), ix_(:)));

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

% if original move was "birth" then reverse move must have children.
assert(~isempty(ix) || code(2) ~= -1); 

if ~isempty(ix)
    transition_prob = 0.5; % chance for deletion
    if (code(2) == 0 && rand < transition_prob) || code(2) == -1
        % erase empty node 
        u = ix(ceil(length(ix)*rand));  % uniformly choose a node to delete 
        transition_prob = transition_prob*1/length(ix);
       
        tree(u,1) = -1; % delete u
        code(2:3) = tree(u,2:3); % record birth/death times for reverse
        code(4) = u;

        was_dead = (tree(u,3) < F); % was deleted node dead?

        ll_diff = 0;
        
        % it's more likely to have less live cells
        if ~was_dead, 
            % record the reads that belong to u
            mv.reads_in_u = find(t==u);
            mv.child_seq = sequences(u,:);
            
            % Associate reads of u to v.        
            t(mv.reads_in_u) = v;
        end
        return;
    end
    % decided to generate a child
    transition_prob = 1-transition_prob;    
end

%   Prepare a proposal to add an empty child:

ll_diff = 0;
% transition_prob = transition_prob;

% draw birth time
birth_time = code(2);
constraint = tree(v,3)-tree(v,2);
if birth_time == 0
    birth_time = tree(v,2) + mod(exprnd(1/birth_rate), constraint); % conditional probability
end

% draw death time (given birth time)
death_time = code(3);
if death_time == 0
    death_time = birth_time + exprnd(1/death_rate);
end
death_time = min(death_time,F); 

% draw sequence of child
B = size(Q,1); L = size(sequences,2);
assert(size(Q,2) == B*L);
parent_ix = double(sequences(v,:)) + (0:B:(B*L-B));
child_seq = int16(many_mnrnd(Q(:,parent_ix),double(mv.child_seq),false));

% no need to take the above component into account in the transition probability
% it cancels out with a likelihood term.  
%%%

%TODO: allocate in advance space for the tree - it takes a lot of time to append.
tree(end+1,:) = [v birth_time death_time];


ll_diff = ll_diff - (birth_rate+death_rate)*(death_time-birth_time);

transition_prob = transition_prob ...
    *exp(-birth_rate*(birth_time-tree(v,2)))/(1-exp(-birth_rate*constraint)) ....
    *exp(-death_rate*(death_time-birth_time));

code(2) = -1;  % mark a "birth" event

% added to support live cells:
 if death_time ==  F,     
    % read assignment likelihood: it's more likely to have less live cells
    ll_diff = ll_diff + length(t)*log(nAlive/(nAlive+1));
    code(3) = 1;
    
    % if parent has reads, move some of them to child
    reads_in_v = find(t==v);

    if isempty(reads_in_v), assert(isempty(mv.reads_in_u)); return; end
    
    seqs = codons2seqs([sequences(v,:); child_seq]);
    mut = (seqs(1,:) ~= seqs(2,:));
    logR = log(R);
    LLs = zeros(2, length(reads_in_v));
    if sum(mut) > 0
        for i=1:length(reads_in_v)
            for j=1:2
                LLs(j,i) = sum(logR(mysub2ind(4, seqs(j,mut), reads(reads_in_v(i),mut),1)));
            end
        end
    end
    % set the forced draw
    forced = zeros(1,length(t));
    if code(3) ~= 0 % this is a reverse move
        forced(reads_in_v) = 1;
        forced(mv.reads_in_u) = 2;
    end
    forced = forced(reads_in_v);
    [t_ prob] = many_lmnrnd(LLs, forced);
    t(reads_in_v(t_==2)) = size(tree,1);
    code(4) = ~isempty(find(t_==2, 1));
        
    % change ll_diff and transition_prob accordingly.
    ll_diff = ll_diff + sum([-1 1]*LLs(:,t_==2));
    transition_prob = transition_prob * prod(prob);    
    
 else 
    code(3) = 0;
 end

end


function sanity(tree)
    
    for u=2:size(tree,1);
        if tree(u,1) == -1, continue; end
        assert(tree(tree(u,1),1) ~= -1);
        assert(tree(u,2)>tree(tree(u,1),2));
        assert(tree(u,2)<=tree(tree(u,1),3));
        assert(tree(u,2)<=tree(u,3));
    end
end