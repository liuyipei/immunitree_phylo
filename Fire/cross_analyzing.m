function [B D] = cross_analyzing(silent, c, save_figures)
%% distance between silent mutation graphs
%close all
if nargin < 3, save_figures = false; end
classes = {'variants_germ', 'clones_germ', 'clones_tree'};
figure;
silent = silent.(c);
%load('~/scratch/data/lymph/t5.fa_results/silent_vectors'); % struct: silent
X = [silent.(classes{1})(:,3:5) silent.(classes{2})(:,3:5) silent.(classes{3})(:,3:5)];
%X(:, [1 2 8 14 18]) = [];
D = pdist(X', 'cityblock');
D = squareform(D);
[a b] = meshgrid(1:size(D,1));
%D(a<b) = 0;
% iy = [1 6 11 15];
iy = [1 4 7 10];
for i=1:3
    classes{i}(classes{i} == '_') = '-';
    for j=1:3
        A = D( iy(i):iy(i+1)-1, iy(j):iy(j+1)-1 );
        B(i,j) = sum(A(:));
        if i==j
            B(i,j) = B(i,j)/(0.5*(length(A(:))- size(A,1)));
        else
            B(i,j) = B(i,j)/length(A(:));
        end

    end
end
clf%figure;
imagesc(D);%, [0 0.5]);
plot_class_lines(iy, 0, 'm3', true)
set(gca, 'XTickLabel', classes);
 set(gca, 'XAxisLocation', 'top');
plot_class_lines(iy, 1, 'm3', true)
%set(gca, 'YTickLabel', classes);
% TODO: I don't understand why I can't change font size
text(0.2*ones(1,9), 1:9, ('smmsmmsmm')', 'fontsize', 13, 'FontWeight', 'bold');
h = text(1:9, 9.5*ones(1,9), ...
    ({'spl', 'mes', 'med', 'spl', 'mes', 'med', 'spl', 'mes', 'med'}),...
    'fontsize', 13, 'FontWeight', 'bold');
set(h, 'rotation', 90);
set(h, 'HorizontalAlignment', 'right');

h = text(0*ones(1,3), [2 5 8], classes, 'fontsize', 13, 'FontWeight', 'bold','rotation', 90, 'HorizontalAlignment', 'center');
set(gca, 'YTickLabel', '');
set(gca, 'fontsize', 13);
set(gca, 'FontWeight', 'bold');

colorbar
title(sprintf('L1 distance between %c-region silent/non-silent hotspots graphs',c));
if save_figures
    saveas(gcf, sprintf('~/www/research/donors/L1_distance_between_%c_region_hotspots_graph.jpg', c));
    saveas(gcf, sprintf('~/www/research/donors/fig/L1_distance_between_%c_region_hotspots_graph.fig', c));
    fprintf('Figures Saved.\n');
end
end


% V =
% 
%     0.1050    0.2348    0.1253
%     0.2348    0.1306    0.2085
%     0.1253    0.2085    0.0439


% J =
% 
%     0.1812    0.3256    0.1501
%     0.3256    0.1924    0.2875
%     0.1501    0.2875    0.0747


%%

function test()
%%
close all
load('~/scratch/data/lymph/t5.fa_results/silent_vectors');
cross_analyzing(silent, 'V', true)
cross_analyzing(silent, 'J', true)


%%

load('~/scratch/data/lymph/t5.fa_results/silent_vectors');
classes = {'clones_germ', 'clones_tree', 'variants_germ'};
c = 'V';
silent = silent.(c);
X = [silent.(classes{1}) silent.(classes{2})];
x = X(:,3);
y = X(:,9);

%%
x= reshape(x,2,[]);
y= reshape(y,2,[]);
figure;
subplot(2,1,1);
bar(x(1,:)-y(1,:), 'g')
    xlim([48 108])
ylim([-0.005 0.015]);
title('silent germ-tree');

subplot(2,1,2);
bar(x(2,:)-y(2,:), 'y')
    xlim([48 108])
ylim([-0.005 0.015]);
title('non silent germ-tree');
%%
figure
l = 0;
for alpha=[0 0.55:0.05:0.95]
%for alpha=1:1:10
    l = l+1;
    subplot(5,2,l);
%    z = (1+alpha)*y - alpha*x;
    z = (y - alpha*x);%/(1-alpha);
    z = reshape(z,2,[]);

    %     x = [sum(reshape([clone_mutations(ix).J_germ_silentspots], NJ, [])'); 
%         sum(reshape([clone_mutations(ix).J_germ_nonsilentspots], NJ, [])')];
    
%    x = x /sum(x(:)); %./ repmat(sum(x,2), 1, size(x,2));
%    z(s,:) = x(:);
    
%    bar(sum(z),'g');
    bar(z(1,:),'g');
    hold on
    bar(-z(2,:), 'y');
    xlim([48 108])
%    ylim([-0.1 0.1]);
    title(sprintf('alpha = %.2f  norm(z)=%.4f', alpha, norm(z,2)));
    
end

end



