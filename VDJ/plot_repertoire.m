function plot_repertoire(rep, V, D, J)
    C = size(rep,2)*ones(1,size(rep, 3));
    rep = reshape(rep, size(rep, 1), []);
    imagesc(rep);
    colorbar;
    plot_class_lines(C);
    ylabel('V');
    xlabel('J,D');
end



function plot_class_lines(C, horizontal)
    if nargin<2, horizontal = 0; end
    % C is a sorted array for labels
    hold on
    chr_ = cumsum(C)+1;
    chr = [1 chr_(1:end-1)];
    chr_middle = ceil((chr+chr_)/2);

    if horizontal == 0        
        for i=1:length(chr)
            plot([chr(i) chr(i)], ylim, 'm-', 'linewidth', 3); 
        end
        set(gca, 'XTick', chr_middle);
        set(gca, 'XTickLabel', 1:length(chr));
    end
    if horizontal == 1        
        for i=1:length(chr)
            plot(xlim,[chr(i) chr(i)], 'm-', 'linewidth', 2); 
        end
        set(gca, 'YTick', chr_middle);
        set(gca, 'YTickLabel', 1:length(chr));
    end
    
    hold off
end