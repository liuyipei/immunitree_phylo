function [X ix depth] = order_tree(T, ix)
    if nargin<2, ix = 1:size(T,1); end
    % get ancsestors table
    X = zeros(length(ix), 1);
    depth = zeros(length(ix), 1);
    for i=1:size(X,1)
        x = ancsestors(ix(i),T(:,1));
        X(i,1:length(x)) = x;
        depth(i) = length(x);
    end
    [X ord]= sortrows(X);
    ix = ix(ord);
    depth = depth(ord);
end


function x = ancsestors(i, T)
    x = [];        
    while(i>0)
        x = [i x];
        i = T(i,1);
    end
    if i == -1, x = []; end
end