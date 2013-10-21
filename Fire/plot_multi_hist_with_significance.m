function plot_multi_hist_with_significance(X, labels, tl, intervals)

multi_hist_with_significance(X, labels, tl, true, intervals, 1, 'random_effects', 'random_effects');
[P H] = multi_hist_with_significance(X, labels, tl, true, intervals, 1, 'ttest', 'mytest');

% this part prints text on the picture.  
for j=1:length(P)
    subplot(length(P), labels.num+1, j*(labels.num+1));
    print_text_on_image(P{j});
end

subplot(length(P), labels.num+1, 1:labels.num);
legend(labels.names, 'location', 'SouthEast');

end