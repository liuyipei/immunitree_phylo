function [W a] = absorb_leaves(W, a)

    [~, ~, ~, src_, rep_] =  play_with_clone(W);
    src = sub2ind([4 6], src_, rep_); % source + replicate information


    nNodes = length(a.tree);
    lonely_leaves = find(a.tree(:,3) == 1);
    lonely_reads = find(ismember(a.t, lonely_leaves));
    read_del = []; 
    node_del = [];
    for i=lonely_reads' % 1:length(lonley_reads)
        node = a.t(i);
        parent_node = a.tree(node, 1);
        reads_at_parent_node = (a.t == parent_node);
        if ~isempty(find( src(reads_at_parent_node) == src(i) , 1))
            % read needs to be removed!
            read_del = [read_del i];
            node_del = [node_del a.t(i)];
        end
    end
    W(1+read_del) = [];
    a.t(read_del) = [];
    assert(length(a.t)+1==length(W));
    a.tree(node_del,1) = -1;
    [a.tree, a.sequences, a.t] = clean_tree(a.tree, a.sequences, a.t);
    for i=1:length(a.t)
        W(i+1).Header(end-2:end) = sprintf('%03d', a.t(i));
    end
%     their_parents = T(lonely_leaves,1);
%     lonely_replicase = replicas(lonely_reads);

end
