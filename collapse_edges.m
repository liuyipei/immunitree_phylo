function b = collapse_edges(a)
    
    parent = a.tree(:,1);
    sequences = a.sequences;
    count = hist(a.t,1:size(a.tree,1))'; 
    t = a.t;

    for k=find(count == 0)'
        if parent(k) ~= 0
            child = find(parent==k);
            if length(child) == 1
                % erase k
                parent(child) = parent(k);
                a.tree2(parent(k),3) = max(a.tree2(parent(k),3), a.tree2(k,3)); 
                parent(k) = -1;  
            end
        end
    end
    a.tree2(:,1) = parent;
    parent = [parent count];
    [parent sequences t] = clean_tree(parent, sequences, t);
    parent = fill_counts(parent);

    b.tree = parent;
    b.sequences = sequences;
    b.t = t;
    b.tree2 = clean_tree(a.tree2);
end

