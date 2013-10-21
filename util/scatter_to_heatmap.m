function z = scatter_to_heatmap(x,y, limx, limy, normalize, clrbar)
    if nargin<5 clrbar = false; end
    nxbins = length(limx)-1;
    nybins = length(limy)-1;
%     limx = 0:(max(x)/nxbins):max(x);
%     limy = 0:(max(y)/nxbins):max(y);
    [~, binx] = histc(x, limx);
    [~, biny] = histc(y, limy);
    binx = min(binx,nxbins);
    biny = min(biny,nybins);
    iz = sub2ind([nxbins,nybins], binx, biny);
    z = hist(iz, 1:nxbins*nybins);
    color_map = [0 1.03.^(1:255)];
    if normalize
        color_map = color_map/color_map(end);
        z = z/sum(z);
    end
    [~,z] = histc(z, color_map);
    
    z = reshape(z, nxbins, nybins);
    z = flipud(z');
%     [cs h] = contour(limx(2:end), limy(2:end), z, 2.^(2:10));
%     clabel(cs, h, [16, 64, 256, 1024]);
    imagesc(z, [1 256]);

    % ticks...  x and y are switched, and x is up side down
    ticks = [1 5 10];
    xticks = (limx(1:end-1)+limx(2:end))/2;
    yticks = (limy(1:end-1)+limy(2:end))/2;
    xlabels = cell(length(ticks),1);
    ylabels = cell(length(ticks),1);
    for t=1:length(ticks) 
        if xticks(ticks(t)) < 10
            xlabels{t} = sprintf('%.1f', xticks(ticks(t)));
        else
            xlabels{t} = sprintf('%.0f', xticks(ticks(t)));
        end
        if yticks(ticks(end-t+1)) < 10
            ylabels{t} = sprintf('%.1f', yticks(ticks(end-t+1)));
        else
            ylabels{t} = sprintf('%.0f', yticks(ticks(end-t+1)));            
        end        
    end   
    set(gca, 'XTick', ticks);
    set(gca, 'XTickLabel', xlabels);
    set(gca, 'YTick', ticks);
    set(gca, 'YTickLabel', ylabels);

    if clrbar
        h = colorbar;
        yticks = get(h, 'YTick');
        for t=1:length(yticks)
            color_map_labels{t} = sprintf('%1.1g', color_map(yticks(t)));
        end
        %    colorbar('YTickLabel', num2cell(color_map([1 64 128 192 256])));
        colorbar('YTickLabel', color_map_labels);
    end
end

