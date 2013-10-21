function st = order_nodes(st)
    T = st.tree(:,1);
    M = length(T);
    DG = sparse(1+T, 2:M+1, true, M+1, M+1);
    if M == 1, 
        order = 1;
    else
        order = graphtopoorder(DG(2:end, 2:end));
    end
    map = zeros(M,1);
    map(order) = 1:M;
    st.t = map(st.t);
    st.tree = st.tree(order,:);
    st.sequences = st.sequences(order,:);
    ix = st.tree(:,1) > 0;
    st.tree(ix,1) = map(st.tree(ix,1));
end
