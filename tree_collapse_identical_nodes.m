% For every node, try to merge it into another with the same sequence,
% ignoring small differences in gapping differences that do not affect the
% actual represetned sequence. Each node attempts to join the member of its
% equal-sequence equivalence class which is most ancestral in some fixed
% order which respects the topology of the tree. Finally, we assert the
% restriction that the destination node must have died no earlier than the
% selected source node. All sequences are transferred from the source to
% the destination. Finally, the children nodes of the destroyed source node
% are NOT given to the destination (consider for example, what if the
% destination node is not a strict ancestor, and was in fact born later?)
% The children of the destroyed source node are instead given to the
% parent of the destroyed source node.

function [tree_, sequences_, t_, reads_, collapse_status] = ...
    tree_collapse_identical_nodes(tree_, sequences_, t_, reads_, F_)

T = tree_(:,1);   
M = length(T);

% Phase one - compute unnormalized conditionals:
% go over tree in reverse topological ordering
% [~,order] = sort(tree_(:,2));

evidence = nodes_with_evidence(tree_,t_);

non5_sequences_arr = zeros(size(sequences_));
for i = 1:size(non5_sequences_arr,1)
    curr_non5_seq = non5(sequences_(i,:));
    non5_sequences_arr(i,1:length(curr_non5_seq)) = curr_non5_seq;
end
[~, raw_from_uniq, uniq_from_raw] = unique(non5_sequences_arr, 'rows'); % these are node sequencse

for i = 1:size(non5_sequences_arr, 1)
   
    x = uniq_from_raw(i); % x is the unique-class
    j = raw_from_uniq(x); % j is the old rep of the class
    if i == j
        continue
    end
    
    if tree_(j,2) > tree_(i,2) || (i == 1 && j ~= 1)
        % all final representatives will have entered this if-check
        raw_from_uniq(x) = i; % i is the new representative -- the earliest born
    else % j was born earlier
    end
    % also update death time of the current best representative
    tree_(raw_from_uniq(x),3) = max(tree_([i j],3));   
        % intermediate-reps' death times 
        % are modified by accident,
        % but they are getting deleted
        % at the end any ways
end


collapse_status = struct(...
    'nodes_merged', 0, ...
    'children_moved', 0, ...
    'reads_moved', 0,...
    'readgaps_realigned', 0);

t_ = raw_from_uniq(uniq_from_raw(t_)); % move reads of the to-be-deceased

non_root = 2:size(tree_); % move the children of the to-be-deceased
tree_(non_root, 1) = raw_from_uniq(uniq_from_raw(tree_(non_root,1)));

representative = (raw_from_uniq(uniq_from_raw) == (1:length(uniq_from_raw))'); % move the children of the deceased
tree_(~representative,1) = -1; % finally delete the to-be-deceased

assert(isempty(intersect(find(~representative), T(1,:))))
assert(sum(T(:,1) == 0)==1)
assert(all(T(T(:,1)>0,1) >= 0)); % that the parent of any non-deleted node is also not deleted.
assert(all(tree_(unique(t_),3)==F_));
[tree_, sequences_, t_] = clean_tree(tree_, sequences_, t_); % clears the deceased nodes
assert(all(tree_(unique(t_),3)==F_));
end

function non5= non5(x)
    non5= x(x~=5);
end