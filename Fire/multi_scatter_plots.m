function multi_scatter_plots(X, pairs, labels, nBins, flds, same, normalize)
if ~exist('normalize', 'var'), normalize = false; end

figure;
f = 0;
for i=1:size(pairs,1)
    ax = zeros(1,labels.num);
    y = X(:,pairs(i,2)); x = X(:,pairs(i,1));
    if same
        xymax = max(max(x),max(y));
        limx = 0:(xymax/nBins):xymax;
        limy = limx;     
    else
        nxbins = nBins;
        nybins = nBins;
        limx = 0:(max(x)/nxbins):max(x);
        limy = 0:(max(y)/nybins):max(y);
    end
    
    iy = ~isnan(y) & ~isnan(x);       
    for s=1:labels.num
        f = f+1;
        ix = (labels.l == s)  & iy;
        ax(s) = subplot(size(pairs,1), labels.num, f);        
        z = scatter_to_heatmap(x(ix), y(ix), limx, limy, normalize, f==1);
        if i == size(pairs,1), xlabel(labels.names{s}); end
        if s == 3, title([flds{pairs(i,2)} ' vs. ' flds{pairs(i,1)}] ); end
    end
    linkaxes(ax, 'xy');
end



end