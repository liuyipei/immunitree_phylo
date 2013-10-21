function T = fill_counts(T)
    M = size(T,1);
    T3 = T(:,2);
    % go over tree in reverse topological ordering
    DG = sparse(1+T(:,1), 2:M+1, true, M+1, M+1);
    order = graphtopoorder(DG(2:end, 2:end));
    rev_order = fliplr(order);
    for k=rev_order(1:end-1)
        T3(T(k,1) ) = T3( T(k,1) ) + T3(k);
    end
    T = [T T3];
end