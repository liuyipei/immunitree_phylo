function [T, sequences, t] = clean_tree(T, sequences, t)
    % renumber the nodes that are not deleted
    % maintain same node order as before.
    T_ = T; 
    ix = T(:,1) >= 0;
    M = sum(ix);
    map = zeros(size(T,1),1);
    map(ix) = 1:M;
    T = T(ix,:);
    if nargin >= 2
        assert(nargout >= 2);
        sequences = sequences(ix,:);
    end
    ix = T(:,1) > 0;
    T(ix, 1) = map(T(ix,1));
    if nargin>2 && ~isempty(t), t = map(t); end
end

