function analyze_trees()
%%  Loading prerequisites
clear
cd ~/JVL/src/phylo/Fire
addpath('.');
addpath('..');
addpath('../VDJ/');
addpath('../util/');

%%  Starting from a .mat file.  
%clear;
set_id = 'variants';
if strcmp(set_id, 'variants'), set_dir = 'nonclones'; else set_dir = set_id; end
with_tree = strcmp(set_id, 'clones');

tree_dir = sprintf('/afs/cs/u/joni/scratch/data/lymph/t5.fa_results/%s/', set_dir)
%tree_dir = sprintf('/afs/cs/u/joni/scratch/data/lymph/t5.fa_results_official/%s/', set_dir)
%tree_dir = sprintf('/afs/cs/u/joni/scratch/data/lymph/t5.fa_results_pooled/%s/', set_dir)
load([tree_dir 'all_clones.mat']);
labels = label_reads('source', tree_dir);
if strcmp(set_id, 'variants'), labels.num = 4; end
output_dir = sprintf('~/www/research/donors/%s/', set_dir)
save_figures = false
% load([tree_dir 'fastas2.mat'])
%%  Figure 1b
% calculate the overlap between the different sources for each patient
[codes singles doubles triplets] = get_tuple_codes();

patient_stats = zeros(10,1);
overlap_stats = zeros(10,16);
% For each patient
for j=1:10
    ix = find(clone_ids(:,1) == j);
    
    % how many clones
    nClones = length(ix);
    
    % distributions over the lengths
    %medLength = median(clone_ids(ix,4));
    
    src = clone_src_(ix,:) > 0; % change this to 1 to consider clones with at least 2 reads from each source
    src = src(:,1)*8 + src(:,2)*4 + src(:,3)*2 + src(:,4);
    h = hist(src, 0:15);   
   
%    patient_stats(j,:) = [nClones medLength];
    patient_stats(j) = nClones;
    overlap_stats(j,:) = h;
end

labels_ = [];
for i=1:length(doubles)
    ix = find(codes(doubles(i),:)); 
    labels_ = [labels_; {sprintf('%d%d', ix(1), ix(2))}];
end

% show figure showing all clone sharing among sources
%close all
figure;
pos = get(gcf, 'position');
pos(3:4) = [700 1300];
set(gcf, 'position', pos);
ind = 1;
L = 10;
if strcmp(set_id, 'variants'), rows = L/2; cols = 2; else cols = 4; rows = L;  end

for j=1:L
    subplot(rows,cols,ind);
    h = bar([overlap_stats(j,singles); 0 0 0 0]);    
    color_the_bars(h);
    xlim([0.5 1.5]);
    junk = ylim;
    ylim([-0.0001 junk(2)]);
    set(gca, 'xtick', []);
    box(gca, 'on')
    ylabel('No. clones');
    
    if j == 1, 
        title('pure clones', 'fontsize', 15); 
    end
    if j == L
        legend(labels.names, 'location', 'SouthOutside'); 
    end    
    
%    title(sprintf('patient %d, total %d %s', j, patient_stats(j,1), set_id));
    ind = ind+1;

    if strcmp(set_id, 'clones')
        subplot(rows,cols,ind);
    %    imagesc(reshape(overlap_stats(j,doubles)/patient_stats(j,1), 4, 4), [0 1]);
        h = bar([overlap_stats(j,doubles); zeros(1,length(doubles))]);        
        color_the_bars(h);
        xlim([0.5 1.5]);
        junk = ylim;
        ylim([-0.01 junk(2)]);
        set(gca, 'xtick', []);
        box(gca, 'on')

        if 0
            xlim([0 7]);
            set(gca,'XTick',1:6)
            set(gca,'XTickLabel',labels_)
        end
        if j==1
            title('shared by 2', 'fontsize', 15); 
        end
        if j == L
            legend('blood-spleen', 'blood-mes', 'blood-med', 'spleen-mes', 'spleen-med', 'mes-med', 'location', 'SouthOutside');
        end

        ind = ind+1;

        subplot(rows,cols,ind);
        h = bar([overlap_stats(j,triplets); 0 0 0 0 ]);        
        color_the_bars(h);
        xlim([0.5 1.5]);
        junk = ylim;
        ylim([-0.0001 junk(2)]);
        set(gca, 'xtick', []);
        box(gca, 'on')

        if j == 1, 
            title('shared by 3', 'fontsize',15); 
        end
        if j == L
            legend({'all but blood', 'all but spleen', 'all but mes', 'all but med'}, 'location', 'SouthOutside');
        end
        ind = ind+1;
        
        if cols == 4
            subplot(rows,cols,ind);
            vn1 = [sum(overlap_stats(j, codes(:,2))) sum(overlap_stats(j, codes(:,3)))  sum(overlap_stats(j, codes(:,4))) ];
            vn2 = [sum(overlap_stats(j, codes(:,2)&codes(:,3)),2) ...
                       sum(overlap_stats(j, codes(:,2)&codes(:,4)),2) ...
                       sum(overlap_stats(j, codes(:,3)&codes(:,4)),2) ...
                       sum(overlap_stats(j, codes(:,2)&codes(:,3)&codes(:,4)),2)];
        %    vn2 = [overlap_stats(j, doubles(4:6))+overlap_stats(j, triplets(1))  overlap_stats(j, triplets(1))];
            if vn1(1) == 0, vn1(1) = vn1(1)+1; end
            venn(vn1, vn2, 'FaceColor',{'y','c', 'b'});%, 'ErrMinMode', 'ChowRodgers');     
            axis equal
            axis off
            if j == L, title('Venn', 'fontsize', 15);  
                legend({'spleen', 'mes', 'med'}); 
            end
            ind = ind+1;
        end
    end
end

if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'clone_sharing.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'clone_sharing.fig'));
    fprintf('Figures saved.\n');
end



%% Figure 1a
% Alternative plot for clone counts -- Variants
L = 10;
cols = L; 
if strcmp(set_id, 'variants')
    rows = 1;
    height = 200;
else
    rows = 3;
    height = 670;    
end
figure;
subplot(rows,1,1)
h = bar(overlap_stats(:,singles));
%title('pure clones')
legend(labels.names(1:4), 'location', 'EastOutside', 'Fontsize', 14);
%xlabel('Donors');
%ylabel(sprintf('No. of %s', set_id));
color_the_bars(h);
set(gcf, 'position', [321 1286 1600 height]);
xlim([0.5 10.5]);
set(gca, 'xtick', []);
if save_figures
    filename = sprintf('%s_counts', set_id);
    saveas(gcf, sprintf('%s/%s.jpg', output_dir, filename));
    saveas(gcf, sprintf('%s/fig/%s.fig', output_dir, filename));
    fprintf('Figures saved.\n');
end

%%
%  Figure 1b - alternative
% show figure showing all clone sharing among sources
%close all, figure;

for j=1:L
    ind = cols+j;
    if strcmp(set_id, 'clones')
        subplot(rows,cols,ind);
    %    imagesc(reshape(overlap_stats(j,doubles)/patient_stats(j,1), 4, 4), [0 1]);
        h = bar([overlap_stats(j,doubles); zeros(1,length(doubles))]);        
        color_the_bars(h);
        xlim([0.5 1.5]);
        junk = ylim;
        ylim([-0.01 junk(2)]);
        set(gca, 'xtick', []);
        box(gca, 'on')
        title(sprintf('Donor %d', j), 'fontsize', 15);
        
         if j == L
             legend('blood-spleen', 'blood-mes', 'blood-med', 'spleen-mes', 'spleen-med', 'mes-med', 'location', 'EastOutside', 'Fontsize', 14);
         end

        ind = ind+cols;

        subplot(rows,cols,ind);
        h = bar([overlap_stats(j,triplets); 0 0 0 0 ]);        
        color_the_bars(h);
        xlim([0.5 1.5]);
        junk = ylim;
        ylim([-0.0001 junk(2)]);
        set(gca, 'xtick', []);
        box(gca, 'on')
        title(sprintf('Donor %d', j), 'fontsize', 15);

         if j == L
             legend({'all but blood', 'all but spleen', 'all but mes', 'all but med'}, 'location', 'EastOutside', 'Fontsize', 14);
         end
        ind = ind+cols;
        
        if rows== 4
            subplot(rows,cols,ind);
            vn1 = [sum(overlap_stats(j, codes(:,2))) sum(overlap_stats(j, codes(:,3)))  sum(overlap_stats(j, codes(:,4))) ];
            vn2 = [sum(overlap_stats(j, codes(:,2)&codes(:,3)),2) ...
                       sum(overlap_stats(j, codes(:,2)&codes(:,4)),2) ...
                       sum(overlap_stats(j, codes(:,3)&codes(:,4)),2) ...
                       sum(overlap_stats(j, codes(:,2)&codes(:,3)&codes(:,4)),2)];
        %    vn2 = [overlap_stats(j, doubles(4:6))+overlap_stats(j, triplets(1))  overlap_stats(j, triplets(1))];
            if vn1(1) == 0, vn1(1) = vn1(1)+1; end
            venn(vn1, vn2, 'FaceColor',{'y','c', 'b'});%, 'ErrMinMode', 'ChowRodgers');     
            axis equal
            axis off
            if j == L, 
                legend({'spleen', 'mes', 'med'}, 'Fontsize', 14, 'location', 'EastOutside'); 
            end
            ind = ind+1;
        end
    end
end

if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'clone_sharing.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'clone_sharing.fig'));
    fprintf('Figures saved.\n');
end


%% Alternative plot for clone counts -- Sharing for clones

figure;
subplot(3,1,2);
h = bar(overlap_stats(:,doubles));
color_the_bars(h);
legend('blood-spleen', 'blood-mes', 'blood-med', 'spleen-mes', 'spleen-med', 'mes-med');
xlabel('donors');
ylabel(sprintf('No. of %s', set_id));
set(gcf, 'position', [321 1286 1239 700]);
title('clones shared by two tissues');

subplot(3,1,3)
h = bar(overlap_stats(:,triplets));
title('clones shared by 3 tissues');
legend(labels.names(1:4));
xlabel('donors');
ylabel(sprintf('no. of %s', set_id));
set(h(1),'FaceColor',[1 0 0], 'DisplayName', 'spleen-mes-med');
set(h(2),'FaceColor',[1 1 0], 'DisplayName', 'blood-mes-med');
set(h(3),'FaceColor',[0 1 1], 'DisplayName', 'blood-spleen-med');
set(h(4),'FaceColor',[0 0 1], 'DisplayName', 'blood-spleen-mes');

if save_figures
    filename = sprintf('%s_counts', set_id);
    saveas(gcf, sprintf('%s/%s.jpg', output_dir, filename));
    saveas(gcf, sprintf('%s/fig/%s.fig', output_dir, filename));
    fprintf('Figures saved.\n');
end



%% Distance from germline plot -- Clones only!
% For each patient show:
% precentage of reads from each source with dist = X from germline
close all

flds = {'height',       'max mutational height';    ...
        'nNodes',       'No. of nodes';             ...
        'medMut',       'median mutation distance'; ...
        'dist_to_root', 'dist. germline to tree root'; ...
        'maxMut',       'max mutation distance'; ...
        'len',          'germline length (after trimming reads)'; ...
        'avg_height',   'avg. height'; ...
        'len_n1',       'length N1'; ...
        'len_n2',       'length N2'; ...
        'med_depth',    'med. depth over nodes'; ...
        'med_depth_leaves', 'med. mut. depth of leaves';...
        'med_depth_nonempty', 'med. mut. depth of nonempty nodes';...
        'leaves_frac', 'fraction of leaves';...
        'avg_node_distance', 'avg. dist. between nonempty nodes'; ...
        }

% [P H] = kstest_with_labels([clone_mat.(fld{1})], labels, fld{2}, 'per patient');        
% H     = hist_with_labels  ([clone_mat.(fld{1})], labels, fld{2}, 'clones', 1, normalize, 'separate');

T = length(flds);
X = zeros(N,T);
for i=1:T
    fld = flds(i,:);
    X(:,i) = [clone_mat.(fld{1})];
end

%% Figure 3
% labels_ = labels;
% labels_.l(ix) = 100;

metric = [1 10 3 5];% 13] ;
plot_multi_hist_with_significance(X(:,metric), labels, flds(metric,2),[1 1 1 1 0.1]);
pos = get(gcf, 'position');
pos(3:4) = [1024 768];
set(gcf, 'position', pos);
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'tree_features1.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'tree_features1.fig'));
    fprintf('Figures saved.\n');
end

%% heatmaps against number of nodes
pairs = [2*ones(length(metric),1) metric'];
multi_scatter_plots(X, pairs, labels, 10, flds(:,2), false, false); % normalized = false
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'number_of_nodes_vs_tree_features1.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'number_of_nodes_vs_tree_features1.fig'));
    fprintf('Figures saved.\n');
end

multi_scatter_plots(X, pairs, labels, 10, flds(:,2), false, true); % normalized = true
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'fraction_of_nodes_vs_tree_features1.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'fraction_of_nodes_vs_tree_features1.fig'));
    fprintf('Figures saved.\n');
end

%%  Remove "blood" from labels
labels.l = labels.l-1;
labels.num = labels.num-1;
labels.names = labels.names(2:end);
%% Per-Patient-Data
close all
for l=1:length(metric)
    hist_with_labels(X(:,metric(l)), labels, flds{metric(l),2}, 'clones', 1, 1, 'separate per patient');   
    set(gcf, 'position', [35  6 1211 1491]);    
    legend(labels.names, 'location', 'SouthEast');
    if save_figures
%        saveas(gcf, sprintf('%s/per_patient_data/%s_cdf.jpg', output_dir, flds{metric(l),1}));
        saveas(gcf, sprintf('%s/per_patient_data/fig/%s_cdf.fig', output_dir, flds{metric(l),1}));
        fprintf('Figures saved.\n');
    end
end
%% another set of features
metric = [2 4 6 8 9];
plot_multi_hist_with_significance(X(:,metric), labels, flds(metric,2),1);
pos = get(gcf, 'position');
pos(3:4) = [1024 680];
set(gcf, 'position', pos);
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'tree_features2.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'tree_features2.fig'));
    fprintf('Figures saved.\n');
end

%% heatmaps against number of nodes (with second set of features)
pairs = [2*ones(length(metric),1) metric'];
multi_scatter_plots(X, pairs, labels, 10, flds(:,2), false);
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'number_of_nodes_vs_tree_features2.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'number_of_nodes_vs_tree_features2.fig'));
    fprintf('Figures saved.\n');
end

multi_scatter_plots(X, pairs, labels, 10, flds(:,2), false, true); % last true means it's normalized
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'fraction_of_nodes_vs_tree_features2.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'fraction_of_nodes_vs_tree_features2.fig'));
    fprintf('Figures saved.\n');
end
%% Heatmaps plots:  
% 1) median distance from node to leaf vs dist to germline
% 2) avg distance between nonempty nodes vs  dist to germline

pairs = [4 12; 4 14];
multi_scatter_plots(X, pairs, labels, 10, flds(:,2), true, true)
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'root2germ_vs_med_node2root.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'root2germ_vs_med_node2root.fig'));
end



%%  read distances (not tree based) from germline and from each other --  clones only

close all
X = clone_intra_read(:,[5 15]);
tl = {'read  dist2germ', 'read hamming'}
X(isnan(X)) = 0;
plot_multi_hist_with_significance(X, labels, tl,1);
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'clone_features1.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'clone_features1.fig'));
end
%%  scatter plot read distances from germline and from each other -- clones
figure;
% for each source
for s = 1:labels.num
    ix = (labels.l == s);
    subplot(1,labels.num,s);
    plot(X(ix,1), X(ix,2), 'x');    
    xlabel(labels.names{s});
    ylim([0 45]);
    xlim([0 60]);
end
subplot(1,labels.num,3);
title('avg distance to germline vs avg distance between reads');
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'dist2germ_vs_dist2read.jpg'));
    saveas(gcf, sprintf('%s/%s', output_dir, 'dist2germ_vs_dist2read.fig'));
end

%%
figure;
% for each source
f = 0;
for k=1:10
    for s = 1:labels.num
        ix = (labels.l == s) & (labels.patients == k);
        f = f+1;
        subplot(10,labels.num,f);
        plot(X(ix,1), X(ix,2), 'x');    
    %    xlabel(labels.names{s});
        ylim([0 40]);
        xlim([0 50]);
    end
end
subplot(10,labels.num,3);
title('avg distance to germline vs avg distance between reads, per source,patient');
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'dist2germ_vs_dist2read_per_patient.jpg'));
end

%%  Distance to germline of variants -- variants FIGURE 3a
close all
X = clone_intra_read(:,[6 6]);% 5]);
ix = clone_intra_read(:,6) <= 1;
X(ix,2) = nan;
flds = {'%mutation in V', '%mutation in V in non-naive cells'};%, 'distance to germline'};
plot_multi_hist_with_significance(X, labels, flds, [0.01 0.1]);% 1]);
set(gcf, 'position', [35  6 1211 1491]);    
legend(labels.names(1:4), 'location', 'SouthEast');
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'variant_features1.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'variant_features1.fig'));
%    saveas(gcf, sprintf('%s/per_patient_data/%s_cdf.jpg', output_dir, 'dist_to_germ'));
end

%% Per-patient plot of percent mutation in V
X = clone_intra_read(:,6);
hist_with_labels(X, labels, '%mutations in V', 'variants', 0:0.01:30, 1, 'separate per patient');   
legend(labels.names);
if save_figures
    saveas(gcf, sprintf('%s/per_patient/%s_cdf.jpg', output_dir, 'mutations_in_v'));
    saveas(gcf, sprintf('%s/per_patient/fig/%s_cdf.fig', output_dir, 'mutations_in_v'));
end

%% Per-patient plot of distance to germline
X = clone_intra_read(:,5);
hist_with_labels(X, labels, 'distance to germline', 'variants', 0:1:35, 1, 'separate per patient');   
legend(labels.names);
if save_figures
    saveas(gcf, sprintf('%s/per_patient/%s_cdf.jpg', output_dir, 'distance_to_germline'));
    saveas(gcf, sprintf('%s/per_patient/fig/%s_cdf.fig', output_dir, 'distance_to_germline'));
end

%%  Supp Figure 4: germline features - how much was eaten.

close all
tl = {'length', 'length N1', 'length N2', 'eaten from V', 'eaten from D1', 'eaten from D2', 'eaten from J', };   
X = [clone_ids(:,2) cellfun(@length, clone_n1dn2(:,1:2)) clone_stats(:,4:7)];
plot_multi_hist_with_significance(X, labels, tl, 1);
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'clone_features2.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'clone_features2.fig'));
end


%%  Bushy Chainy

labels = label_reads('source', tree_dir);
Y = [ [clone_triplets.chainy]' [clone_triplets.bushy]'];
Y(Y == -1) =NaN;
%X(:,1) = X(:,1)/100;

plot_multi_hist_with_significance(Y, labels, {'chainy', 'bushy'}, [0.1 1]);
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'bushy_chainy.jpg'));
    saveas(gcf, sprintf('%s/%s', output_dir, 'bushy_chainy.fig'));
end


%%
figure;
flds = {'chainy', 'bushy'};
f = 0;
for i=1:size(Y,2)
    ax = zeros(1,labels.num);
    for s=1:labels.num
        f = f+1;
        ax(s) = subplot(size(Y,2), labels.num, f);
        ix = labels.l == s;
        %z = scatter_to_heatmap(X(ix,metric(i)), X(ix,2));
        %imagesc(z, [0 100]);
        %contour(z);
        scatter(X(ix,2), Y(ix,i), '.');
        if i == size(Y,2), xlabel(labels.names{s}); end
        if s == 3, title(flds{i}); end
    end
    linkaxes(ax, 'xy');
    xlim([0 max(X(:,2))]);
    ylim([0 max(Y(:,i))]);
end
saveas(gcf, 'number_of_nodes_vs_bushy_chainy.jpg');



%%  scatter plot bushy-chainy
figure;
% for each source
for s = 1:labels.num
    ix = (labels.l == s);
    subplot(1,5,s);
    plot(X(ix,1), X(ix,2), 'x');    
    xlabel(labels.names{s});
    ylabel('bushy');
    ylim([0 20]);
    xlim([-0.05 1.05]);
end
subplot(1,5,3);
title('bushy(y) vs. chainy(x)');
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'bushy_vs_chainy.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'bushy_vs_chainy.fig'));
end


%%   Estimate transition matrices
ix = false(16,1); 
ix(mysub2ind(4, 1:4, 1:4, 1)) = true;
%NT = [];

rep = load_repertoire('ihmmune_collapsed');
for c='J'
    mutations = [clone_mutations.(c)];
    L = length(mutations(1).germ_hotspots);
    A = reshape([mutations.germ_hotspots], L, [])';
    
    % hist every column to the values 0:16
    h = hist(A, 0:20);
    h = h(2:17,:);    
    h = reshape(h,4,4,[]); % parent x child x location    
    h = permute(h, [2 1 3]); % child x parent x location
    h = reshape(h, 16, []);
    h = h+2;    
    h(ix,:) = h(ix,:) + 98;
    h = reshape(h, 4, []);
    h_ = h ./ repmat(sum(h),4,1);   
    NT.(c) = reshape(h_, 4, 4, []);    
    for k=1:size(NT.(c),3)
        NT.(c)(:,:,k) = real(NT.(c)(:,:,k)^0.04);
    end
    NT.(c)(NT.(c)<0) = eps;
end





%%  V annotation
annot.FR1 = [1 26];
annot.CDR1 = [27 38]; 
annot.FR2 = [39 55];
annot.CDR2 = [56 65];
annot.FR3 = [66 104];
annot.CDR3 = [105 116];

flds = fields(annot);
annot.ref = zeros(1,annot.CDR3(2));
for k = 1:length(flds)
    ix = annot.(flds{k});
    annot.ref(ix(1):ix(2)) = k;
end
plot_class_lines(annot.ref)
set(gca, 'XTickLabel', flds);

%%  Supp Figure 2,3: V/J-region aligned with silent/non-silent
%close all

annot.FR1 = [1 26];
annot.CDR1 = [27 38]; 
annot.FR2 = [39 55];
annot.CDR2 = [56 65];
annot.FR3 = [66 104];
annot.CDR3 = [105 116];
flds = fields(annot);
annot.ref = zeros(1,annot.CDR3(2));
for k = 1:length(flds)
    ix = annot.(flds{k});
    annot.ref(ix(1):ix(2)) = k;
end

per_clone = false;
if per_clone
    f_str = '_per_clone';
    total_ent = set_id;
else
    f_str = '';
    total_ent = 'mutations';
end
%silent = [];

for c='VJ'
    rep = load_repertoire('ihmmune_collapsed');
    mutations = [clone_mutations.(c)];
    L = length(mutations(1).germ_silentspots);

    % A = get_IGH_alignments(rep, c);
    % A = A(clone_ids(:,2), 3:3:end);
    A = reshape([mutations.support], L, [])';
    A_ = sum(A,1);
    thresh = 0; 
    if strcmp(set_id, 'clones'), thresh = 500; else thresh = 8500; end
    nullA = (A_ < thresh);
    A_(:,nullA) = [];
    if c == 'V', 
        annot.ref(nullA) = []; 
        flds = flds(annot.ref(1):end); 
        annot.ref = annot.ref-annot.ref(1)+1; 
    end
    
    figure;
    if with_tree, r=2; else r=1; end

    ax = subplot(r,1,1);
    x = [sum(reshape([mutations.germ_silentspots], L, [])'); 
        sum(reshape([mutations.germ_nonsilentspots], L, [])')];
    x(:,nullA) = [];


    if per_clone
        x = x./ [A_; A_]; % fraction of clones with mutation there
    else
        x = x / sum(x(:));  % distribution over hotspots
    end

    bar(x(2,:),'g');
    hold on
    bar(-x(1,:), 'y');

%    if ~per_clone
        maxx = max(x(:));   
%         plot( maxx*A_/N, '--', 'linewidth',2);
%         plot(-maxx*A_/N, '--', 'linewidth',2);
%    end

    xlabel('aligned codon position in germline');
    %xlim([46 104]);
    title([c '-region hotspots for non/silent mutations relative to germline']);
    ylabel(sprintf('Fraction of %s', total_ent), 'fontsize', 14);
    colormap summer
    legend({'non-silent', 'silent'});
    set(gcf, 'position', [398 1213 1162 284]);
    if c == 'V'
        plot_class_lines(annot.ref)
        plot_class_lines(annot.ref)
        set(gca, 'XTickLabel', flds, 'fontsize', 11);
    else
        set(gca, 'XTick', 0:2:16, 'fontsize', 10);
    end
    silent_vectors = x(:);

    if with_tree
        B = A ;%.* repmat([clone_mat.nNodes]', 1, size(A,2));
        A_ = sum(B,1);
        A_(:,nullA) = [];
        ax1 = subplot(2,1,2);
        x = [sum(reshape([mutations.silentspots], L, [])'); 
            sum(reshape([mutations.nonsilentspots], L, [])')];
        x(:,nullA) = [];
        if per_clone
            x = x./ [A_; A_]; % fraction of clones with mutation there
        else
            x = x / sum(x(:));  % distribution over hotspots
        end

        hold on
        bar(x(2,:),'g');
        hold on
        bar(-x(1,:), 'y');

%        if ~per_clone
%          maxx = max(x(:));   
%             plot( maxx*A_/N, '--', 'linewidth',2);
%             plot(-maxx*A_/N, '--', 'linewidth',2);
%        end

        xlabel('aligned codon position in germline');
        ylabel(sprintf('Fraction of %s', total_ent), 'fontsize', 14);
    %    xlim([46 104]);
        title([c '-region hotspots for non/silent mutations in tree']);
        set(gcf, 'position', [398 1213 1162 484]);
        if c == 'V'
            plot_class_lines(annot.ref)
            plot_class_lines(annot.ref)
            set(gca, 'XTickLabel', flds, 'fontsize', 11);
        else
            set(gca, 'XTick', 0:2:16, 'fontsize', 10);
        end
        
        linkaxes([ax ax1], 'xy');        
        tree_silent_vectors = x(:);
    end
    ylim([-maxx maxx]);

    if save_figures    
        figure_str = sprintf('%c_region_silent_mutation_hotspots%s', c, f_str);
        saveas(gcf, sprintf('%s/%s.jpg', output_dir, figure_str));
        saveas(gcf, sprintf('%s/fig/%s.fig', output_dir, figure_str));
        fprintf(sprintf('Figures saved in %s.\n', output_dir));
    end

    %  V/J-region aligned with silent/non-silent per source
    ax1 = [];
    if with_tree
        figure;
        ax1 = zeros(1,labels.num);
        for s=1:labels.num
            ix = labels.l == s;
            A_ = sum(B(ix,:),1);
            A_(:,nullA) = [];
            ax1(s) = subplot(labels.num,1, s);

            x = [sum(reshape([mutations(ix).silentspots], L, [])'); 
                sum(reshape([mutations(ix).nonsilentspots], L, [])')];

            x(:,nullA) = [];
            if per_clone
                x = x./ [A_; A_]; % fraction of clones with mutation there
            else
                x = x / sum(x(:));  % distribution over hotspots
            end

            bar(x(2,:),'g');
            hold on
            bar(-x(1,:), 'y');


            title(labels.names{s}, 'fontsize', 11);
            tree_silent_vectors = [tree_silent_vectors x(:)];
            if c == 'V', 
                plot_class_lines(annot.ref); 
                set(gca, 'XTickLabel', flds, 'fontsize', 11);
            else
                set(gca, 'XTick', 0:2:16, 'fontsize', 10);
            end
        end
        subplot(labels.num,1,1); 
        %title([c '-region non/silent mutations in tree']);
        legend({'non-silent', 'silent'});
        subplot(labels.num,1,ceil(labels.num/2)); 
        ylabel(sprintf('Fraction of %s', total_ent), 'fontsize', 14);
        colormap summer
        set(gcf, 'position', [398 1213 1162 884]);
        
        if save_figures
            figure_str = sprintf('%c_region_silent_mutation_hotspots_per_source_tree%s', c, f_str);
            saveas(gcf, sprintf('%s/%s.jpg', output_dir, figure_str));
            saveas(gcf, sprintf('%s/fig/%s.fig', output_dir, figure_str));
            fprintf(sprintf('Figures saved in %s.\n', output_dir));
        end
    end

    figure;
    z =[];% zeros(labels.num, 214);
    ax = zeros(1,labels.num);
    for s=1:labels.num
        ix = labels.l == s;
        A_ = sum(A(ix,:),1);
        A_(:,nullA) = [];

        ax(s) = subplot(labels.num,1, s);
        x = [sum(reshape([mutations(ix).germ_silentspots], L, [])'); 
            sum(reshape([mutations(ix).germ_nonsilentspots], L, [])')];

        x(:,nullA) = [];
        if per_clone
            x = x./ [A_; A_]; % fraction of clones with mutation there
        else
            x = x / sum(x(:));  % distribution over hotspots
        end

        z = [z; x(:)'];

        bar(x(2,:),'g');
        hold on
        bar(-x(1,:), 'y');

        maxx = max(x(:));   

        
        title(labels.names{s}, 'fontsize', 11);
        silent_vectors = [silent_vectors x(:)];
        if c == 'V', 
            plot_class_lines(annot.ref); 
            set(gca, 'XTickLabel', flds, 'fontsize', 11);
        else
            set(gca, 'XTick', 0:2:16, 'fontsize', 10);
        end                
    end
    ylim([-maxx maxx]);
    linkaxes([ax ax1], 'xy');
    subplot(labels.num,1,1); 
    %title([c '-region non/silent mutations root compared to germline']);
    legend({'non-silent', 'silent'});
    subplot(labels.num,1,ceil(labels.num/2)); 
    ylabel(sprintf('Fraction of %s', total_ent), 'fontsize', 14);
    colormap summer
    set(gcf, 'position', [398 1213 1162 884]);

    if save_figures
        figure_str = sprintf('%c_region_silent_mutation_hotspots_per_source_germline%s', c, f_str);
        saveas(gcf, sprintf('%s/%s.jpg', output_dir, figure_str));
        saveas(gcf, sprintf('%s/fig/%s.fig', output_dir, figure_str));    
        fprintf(sprintf('Figures saved in %s.\n', output_dir));
    end

    silent.(c).([set_id '_germ']) = silent_vectors;
    if with_tree, silent.(c).clones_tree = tree_silent_vectors; end
end
%%
silent.(c).variants_germ = silent_vectors;


%%   Clone VJ repertoir.  Create a 'data' structure
data = [];
labels_ = labels;
labels_.l(labels.l == 5) = 1000;

data.ihmmune_collapsed = [clone_ids(:,2) zeros(N,1), clone_ids(:,3)];
data.subject_num = (labels_.patients-1)*4 + labels_.l;
data.trimmed_VJ = zeros(length(labels.l), 4);
if strcmp(set_id, 'clones')
    data.dist2germ = [clone_mat.dist_to_root]';
else 
    data.dist2germ = clone_intra_read(:,5);
end
%%
%[~, ~, big_hist] = collect_VJs(data, get_reads(data, 1:40, []), 'ihmmune_collapsed', 0);

rep = load_repertoire('ihmmune_collapsed');
for k=1:40
    [~, ~, Y] = collect_VJs(data_, get_reads(data_,k), 'ihmmune_collapsed', 0);
    if length(rep.V) > size(Y,1), Y(length(rep.V), length(rep.J)) = 0; end
    Y = Y(:);        
    l(k) = binomial(Y) + sum(Y.*log(X));
end


%%  print trees with colors of replicas %%
% ix should be the indexes of trees from a single source (Say blood)
%ix = find(labels.l == 1)';
ix = find([clone_mat.nNodes]' > 6 & labels.l < 5);
iy = randperm(length(ix));
clones_to_print = 25;
samples_dir = 'samples2';

%rnd = ceil(length(ix)*rand(1,10));
cnt = 0;
if ~exist([tree_dir '/' samples_dir '/'], 'dir')
    system(['mkdir ' tree_dir '/' samples_dir '/']);
end

for i=ix(iy)'
    close all hidden
    if mod(i,1000)==0, fprintf('%d ', i); end
    if sum(clone_src(i,:)>0) == 1, continue; end   % ignore single-replica clones

    [clone G id src replicas st] =  play_with_clone(clone_fasta{i}, states(i));
%    [clone G id src replicas st] =  play_with_clone([tree_dir trees(i).name]);

    cnt = cnt+1;
    h = visualize_tree(st.tree, st.sequences, [st.t replicas], clone, 6, G);
    
    if ~isempty(h)
        filename = trees(i).name;  filename(filename == '_') = ' ';
        h.Nodes(1).Label = [h.Nodes(1).Label sprintf('  %s\nindex: %d ', filename, i) ];
        jpg_file_name = [tree_dir '/' samples_dir '/' trees(i).name(1:end-3) '.2.jpg'];
        saveas(h.hgAxes, jpg_file_name); 
    end  
    show_read_alignment_to_germline(clone, G, replicas, st.t);    
    set(gcf, 'Position', [6 281 1126 246]);
    saveas(gcf, [tree_dir '/' samples_dir '/' trees(i).name(1:end-3) '.1.jpg']);
    if cnt == clones_to_print, break; end
end


%% print clones that have more than one source:
ix = find(sum(clone_src_,2)>10 & sum(clone_src_ > 1, 2) > 1);
iy = 1:45; %randperm(length(ix));
clones_to_print = 25;
samples_dir = 'pure'; % shared

fprintf('Processing %d clones\n', length(ix));

cnt = 0;
if ~exist([tree_dir '/' samples_dir '/'], 'dir')
    system(['mkdir ' tree_dir '/' samples_dir '/']);
end

h = [];

for i=(iy)'%ix'%(iy)'
    i
    close all hidden

    [clone G id src replicas st h] =  play_with_clone(clone_fasta{i}, states(i));

    cnt = cnt+1;

    % TODO: confirm that you are writing the intra_read distance
    % excluding same-replica distances
    str = sprintf('\nread2germ: [%.1f %.1f %.1f %.1f %.1f]', clone_intra_read(i,1:5));
    str = [str sprintf('\nread2read: [%.1f %.1f %.1f %.1f %.1f]', clone_intra_read(i,11:15))];
    
    if ~isempty(h)
        filename = trees(i).name;  filename(filename == '_') = ' ';
%        h.Nodes(1).Label = [h.Nodes(1).Label sprintf(' %s index:%d ', filename, i) str];
        jpg_file_name = [tree_dir '/' samples_dir '/' trees(i).name(1:end-3) '.2.jpg'];
        saveas(h.hgAxes, jpg_file_name); 
    end            

    show_read_alignment_to_germline(clone, G, src, st.t);    
    set(gcf, 'Position', [6 281 1126 246]);
    saveas(gcf, [tree_dir '/' samples_dir '/' trees(i).name(1:end-3) '.1.jpg']);
end


%%  Counts percentage of V that is mutated
X = [clone_mutations.V];
Y = [X.germ_hotspots];
Y = reshape(Y,[],N);
support = sum(Y>0);
Y(Y == 1) = 0;
Y(Y == 6) = 0;
Y(Y == 11) = 0;
Y(Y == 16) = 0;
muts = sum(Y>0);
clone_intra_read(:,1) = (100*muts./support)';



%% spearman correlation between dist to germ and no of nodess.
X = [clone_mat.dist_to_root];
Y = [clone_mat.nNodes];
scatter(X+0.1*randn(size(X)),Y+0.1*randn(size(Y)))
xlabel('dist to germline'); ylabel('no of nodes')
corr(X',Y')
title(sprintf('spearman = %.2f', corr(X', Y', 'type', 'Spearman')));

X = X+1;
h = hist(sub2ind([max(X) max(Y)], X, Y), 1:(max(X)*max(Y)));
imagesc(reshape(h, max(X), max(Y)));
ylabel('distance to germline');
xlabel('number of nodes');
xlim([0.5 15]);
ylim([0.5 30]);
colorbar

ticks = get(gca, 'ytick');
set(gca, 'yticklabel', ticks+1);
title(sprintf('heatmap: no. of clones binned by dist. to germline and no. of nodes\nspearman correlation = %.2f', corr(X', Y', 'type', 'Spearman')));


%%  Obsolete!!!






%%  count number of locations different than germline
% allow one to be different than the consensus

%load([tree_dir 'states.mat'])
load([tree_dir 'fastas.mat'])
%%
% shared = at least two sources with at least two reads each
for i=1:N%ix(1)'    
    if mod(i,1000)==0, fprintf('%d\n', i); end
    [clone, G, read_ix, src_, rep_] =  play_with_clone(clone_fasta{i});       
    [Y total] = view_read_alignment(clone, G.seq);
    if size(clone,1) < 4
        X = sum(total(1:4,:)) >= size(clone,1);
    else
        X = sum(total(1:4,:)) >= size(clone,1)-1;
    end
    clone_intra_read(i,1) = sum(X);
end
%%
ix = (sum(clone_src > 0, 2) > 1);

figure;
subplot(2,1,1)
hist(clone_intra_read(~ix,1), 1:90);
ylabel('no. clones');
title('number of locations consensus different from germline');
xlim([0 60]);

subplot(2,1,2)
hist(clone_intra_read(ix,1), 1:90);
ylabel('no. clones');
title('number of locations consensus different from germline, in shared clones');
xlim([0 60]);



%% multi source clones

% how many *replicates* produced >=2 reads for this clone, for each sample
clone_src__ = zeros(size(clone_src_));
for i=1:N, 
    x = reshape(clone_src(i,:) > 0, [40 6]);
    clone_src__(i,:) = sum(x');
end

% multi-source

just_shared = (sum(clone_src__ > 0, 2)>1);
fprintf('shared between two sources (possibly same patient) = %d\n', sum(just_shared));

shared_no_leakage = (sum(clone_src__ > 1, 2)>1);
fprintf('shared between sources, no leakage = %d\n', sum(shared_no_leakage));

shared_possible_leakage = just_shared & ~shared_no_leakage;
fprintf('shared between sources, possible leakage = %d\n', sum(shared_possible_leakage));

close_to_germline = (clone_intra_read(:,1) <=6);
fprintf('close to germline, shared, no leakage = %d\n', sum(close_to_germline & shared_no_leakage));

%% multi-patient

% how many reads came from each patient
clone_patient  = zeros(N,10);  
clone_patient_ = zeros(N,10); % reads from how many replicates came for each patient
for i=1:size(clone_src_,1), % % clone_src_ , clone_src__ are N x 40
    x= reshape(clone_src_(i,:), 4, []); 
    clone_patient(i,:) = sum(x,1); 
    x= reshape(clone_src__(i,:), 4, []); 
    clone_patient_(i,:) = sum(x,1); 
end

shared_between_patients = (sum(clone_patient_ > 0, 2)>1);
fprintf('shared between two patients = %d\n', sum(shared_between_patients));

shared_patient_no_leakage = (sum(clone_patient_ > 1, 2)>1);
fprintf('shared between two patients, no leakage = %d\n', sum(shared_patient_no_leakage));

shared_patient_possible_leakage = shared_between_patients & ~shared_patient_no_leakage;
fprintf('shared between patients, possible leakage = %d\n', sum(shared_patient_possible_leakage));

fprintf('close to germline, patient shared, no leakage = %d\n', sum(close_to_germline & shared_patient_no_leakage));

%junk = [just_shared shared_possible_leakage close_to_germline];


%%  Find average amount of leakage

p34 = clone_patient(:,3)>0 & clone_patient(:,4)>0;
% focus only on 2 member clones
y = sum(clone_patient,2);

sum(shared_between_patients & y==2)/424
sum(shared_between_patients & y==3)/424
sum(shared_between_patients & y==4)/424
sum(shared_between_patients & y==5)/424

sum(y==2)/3863
sum(y==3)/3863
sum(y==4)/3863
sum(y==5)/3863
%%
s = 3;
H = zeros(10);
for i=1:10,
    for j=1:10
        H(i,j) = sum(clone_patient(:,i)==1 & clone_patient(:,j)==s-1 &  (y==s));
    end
    H(i,i) = sum(clone_patient(:,i)==y & (y==s));
end
noleak = sum(diag(H))
leak = (sum(H(:))-sum(diag(H)))
leak/(45*16) * 60 / noleak
%%
H = zeros(10);
for i=1:10,
    for j=1:10
        H(i,j) = sum(clone_patient(:,i)>0 & clone_patient(:,j)>0);
    end
    H(i,i) = sum(clone_patient(:,i)==y);
end


%%   regular clones - I think this can be erased
clone_src__ = zeros(size(clone_src_));
for i=1:N, 
    x = reshape(clone_src(i,:) > 0, [4 6]);
    clone_src__(i,:) = sum(x');
end
shared_no_leakage = (sum(clone_src__ > 1, 2)>1);
sum(shared_no_leakage)
just_shared = (sum(clone_src_ > 0, 2)>1);
sum(just_shared)
shared_possible_leakage = just_shared & ~shared_no_leakage;
sum(shared_possible_leakage)
close_to_germline = (clone_intra_read(:,1) <=6);
sum(close_to_germline)

junk = [just_shared shared_possible_leakage close_to_germline];






%%
samples_dir = 'shared2';

if ~exist([tree_dir '/' samples_dir '/'], 'dir')
    system(['mkdir ' tree_dir '/' samples_dir '/']);
end
for i=ix'
    [clone, G, read_ix, src_, rep_] =  play_with_clone(clone_fasta{i});    
    
    lbl = [ceil(src_/4) mod(src_-1,4)+1  rep_];
    lbl = num2str(lbl);
    lbl = mat2cell(lbl,ones(size(lbl,1),1), size(lbl,2));
    
    show_read_alignment_to_germline(clone, G, mod(src_-1,4)+1,lbl);
    set(gcf, 'Position', [6 281 1126 246]);
    saveas(gcf, [tree_dir '/' samples_dir '/tree_id_' num2str(i) '.1.jpg']);
end




%% plot figure of stats VDJ mutation rate
close all
labels = label_reads('source', tree_dir);
normalize = 1;

% hist_with_labels(clone_stats(:,1), labels, 'Number of reads', 'clones', 1, normalize);

%figure
mut = reshape([clone_mat.mut], [], 3) ./ repmat([clone_mat.nNodes], 3, 1)';
VDJ = 'VDJ';
tl = {};
for i=1:3
    tl{i} = sprintf('total mutations in %c / length(%c) / nSeq', VDJ(i), VDJ(i));
end
multi_hist_with_significance(mut, label_reads('source', tree_dir),          tl, true, 0.01);






%%  Compute "normalized clonality"
% gather statistics for variants
[~, replicate_id] = max(clone_src>0,[],2);
bin_id = sub2ind([24 10], replicate_id, labels.patients);
% h_variants = hist(replicate_id, 1:24);
% h_variants_ = reshape(h_variants, 4, 6);
h_variants = hist(bin_id, 1:240);
h_variants = reshape(h_variants, 4, 6, 10);

%% [load clones]
% gather statistics for clones
h_clones = zeros(4,6,10);
for i=1:10
    h_clones(:,:,i) = reshape(sum(clone_src(labels.patients==i,:)), 4, 6);
end

%%
h_all = h_clones + h_variants;
normalized_clonality = zeros(0,4);

% for each tissue, total cross-replica pairs
for patient=1:10
h = h_all(:,:,patient);

cross_replica_pairs = zeros(1,4);
same_clone_cross_replica_pairs = zeros(1,4);

for s=1:4
    cnt = 0;
%     for patient=1:10
%         h = h_all(:,:,patient);
        for i=1:6
            for j=(i+1):6
                cnt = cnt + h(s,i)*h(s,j);
            end
        end
%    end    
    cross_replica_pairs(s) = cnt;    

    % go over all clones count intra-clone-cross-replica-pairs
   cnt = 0;
%   for i=1:N
   for i=find(labels.patients == patient)'
        X = reshape(clone_src(i,:), 4, 6);
        for j=1:6
            for k=(j+1):6
                cnt = cnt + X(s,k)*X(s,j);
            end
        end                
    end
    same_clone_cross_replica_pairs(s) = cnt;    
end
normalized_clonality = [normalized_clonality; (same_clone_cross_replica_pairs./cross_replica_pairs)];
fprintf('%.3g   \t', normalized_clonality(end,:));
fprintf('\n');
end

%%  Save to a Latex Table 1 (fraction of reads in clones)

    X = reshape(sum(h_clones,2), 4, 10)'./reshape(sum(h_all,2), 4, 10)';
    X = [X; (sum(sum(h_clones,2),3)./sum(sum(h_all,2),3))'];
    columnLabels = labels.names(1:4);
    rowLabels = {'donor 1', 'donor 2', 'donor 3', 'donor 4', 'donor 5', ...
        'donor 6', 'donor 7', 'donor 8', 'donor 9', 'donor 10', 'all'};
    latex_dir = '/afs/cs/u/joni/JVL/src/latex/BD_tree/';
    latex_table = [latex_dir 'reads_in_clones.tex'];
    
    matrix2latex(X, latex_table, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f');
%%  Save to a Latex Table 2 (normalized clonality)

    X = 10000*[normalized_clonality; 0.000164   	0.000119   	5.01e-05   	0.000376];
    columnLabels = labels.names(1:4);
    rowLabels = {'donor 1', 'donor 2', 'donor 3', 'donor 4', 'donor 5', ...
        'donor 6', 'donor 7', 'donor 8', 'donor 9', 'donor 10', 'all'};
    latex_dir = '/afs/cs/u/joni/JVL/src/latex/BD_tree/';
    latex_table = [latex_dir 'normalized_clonality.tex'];
    
    matrix2latex(X, latex_table, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f');

    
    
    %%  Venn diagram combining all individuals -- obselete!
close all

figure;
subplot(1,3,2:3);
vn1 = sum([sum(overlap_stats(:, codes(:,2)),2) sum(overlap_stats(:, codes(:,3)),2) sum(overlap_stats(:, codes(:,4)),2) ]);
vn2 = sum([sum(overlap_stats(:, codes(:,2)&codes(:,3)),2) ...
           sum(overlap_stats(:, codes(:,2)&codes(:,4)),2) ...
           sum(overlap_stats(:, codes(:,3)&codes(:,4)),2) ...
           sum(overlap_stats(:, codes(:,2)&codes(:,3)&codes(:,4)),2)]);

try       
    [~, S] = venn(vn1, vn2, 'FaceColor',{'y','c', 'b'}, 'ErrMinMode', 'ChowRodgers');     
    for i=1:7
        a = text(S.ZoneCentroid(i,1)-5, S.ZoneCentroid(i,2)-2, num2str(S.ZonePop(i)));
        set(a, 'Color', 'k');
    end
end
axis equal
axis off
title(sprintf('shared %s, aggregated', set_id));    
legend({'spleen', 'LN1', 'LN2'});


% Heat map patient/source number of pure read clusters
subplot(1,3,1);
imagesc([overlap_stats(:,singles) sum(overlap_stats,2)-sum(overlap_stats(:,singles),2)]);
ylabel('patients');
xlabel('sources');
title(sprintf('%s per source', set_id));
colormap hot
colorbar

if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'clone_sharing2.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'clone_sharing2.fig'));
    fprintf('Figures saved.\n');
end
