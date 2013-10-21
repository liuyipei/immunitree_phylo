function DG = get_graph_from_tree(T, weights, undirected)
    if ~exist('weights', 'var') || isempty(weights), weights = true; end
    if ~exist('undirected', 'var'), undirected = false; end
    M = size(T,1);
    DG = sparse(1+T(:,1), 2:M+1, weights, M+1, M+1); 
    DG = DG(2:end, 2:end);       
    if undirected
        DG = tril(DG+DG');
    end
end