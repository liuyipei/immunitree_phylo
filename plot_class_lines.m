function plot_class_lines(C, horizontal, color, show_ticks)

    if nargin<4, show_ticks = true; end
    if nargin<3, color(2) = '1'; color(1) = 'k'; end
    if length(color) == 1, color(2) = '1'; end
    lw = double(color(2)-'0');
    if nargin<2, horizontal = 0; end
    C = C(:);
    % C is a sorted array for labels
    hold on
    if any(diff(C) == 0)
        M = max(C);
        chr = zeros(1,M); 
        for i=1:M, chr(i) = find(C(:,1) == i, 1); end
        chr_ = [chr(2:end) length(C)+1];
    else
        assert(C(1) == 1);
        chr = C(1:end-1);
        chr_ = C(2:end);
%         chr_ = cumsum(C)+1;
%         chr = [1 chr_(1:end-1)];
    end
    chr_middle = ((chr+chr_)/2)-0.5;

    if horizontal == 0        
        for i=1:length(chr)
            plot([chr(i)-0.5 chr(i)-0.5], ylim, sprintf('%c-', color(1)), 'linewidth', lw); 
        end
        if show_ticks
            set(gca, 'XTick', chr_middle);
            set(gca, 'XTickLabel', 1:length(chr));
        end
    end
    if horizontal == 1        
        for i=1:length(chr)
            plot(xlim,[chr(i)-0.5 chr(i)-0.5], sprintf('%c-', color(1)), 'linewidth', lw); 
        end
        if show_ticks
            set(gca, 'YTick', chr_middle);
            set(gca, 'YTickLabel', 1:length(chr));
        end
    end
    
    hold off
end