function heatmap_to_scatter(X, rep)
    L = max(abs(X(:)));
    [I J] = ind2sub(size(X), 1:numel(X));
    ix = X(:) > 0;
    scatter(J(ix), I(ix), 150*abs(X(ix))/L, 'bo', 'filled');
    hold on
    ix = X(:) < 0;
    scatter(J(ix), I(ix), 150*abs(X(ix))/L, 'ro', 'filled');
    xlim([0 size(X,2)+1]);
    ylim([0 size(X,1)+1]);
    
    if nargin==2
        Vs = {rep.V.Header};
        Js = {rep.J.Header};
        h = text(1:length(Vs), 2+size(X,1)+1*(-1+mod(1:length(Vs),2)), ...
            cellfun(@(x) x(4:end), Vs, 'uniformoutput', false));
        set(h, 'rotation', 90);
        set(h, 'fontsize', 8);        
    end
    
end