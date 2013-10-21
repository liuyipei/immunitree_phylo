clear
cd ~/JVL/src/phylo/Fire
addpath('.');
addpath('../../phylo/');
addpath('../../DPtrees/');
addpath('../../VDJ/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));
%% mean distance between reads on each tree
% TODO: silent/nonsilent mutations (problem - no frameshift info! )
clear
%tree_dir = '/afs/cs/u/joni/scratch/data/lymph/t5.fa_results_take12_clones_greater_than_10_reads/';
tree_dir = '/afs/cs/u/joni/scratch/data/lymph/t5.fa_results/';
trees = dir([tree_dir '*.fa']);
N = length(trees);
%%
assert(false);
clone_dist = zeros(N,5);
for i=1:N
    if mod(i,1000)==0, fprintf('%d\n', i); end
    
    % load clone .mat file
    mat_file_name = [tree_dir  trees(i).name(1:end-3) '.mat'];
    Z = load(mat_file_name);    
    [nMuts mut] = annotate_mutations_on_tree(Z.a.tree(:,1), Z.a.sequences);
    T = Z.a.tree;

    % load clone .fa file
    W = fastaread([tree_dir trees(i).name]);    
    W_ = {W(2:end).Header};
    W_ = cell2mat(cellfun(@(x) [x '_'], W_, 'uniformoutput', false));
    [~, ~, src t] = strread(W_,'%s%d%*s%d%*s%d','delimiter','_');
    src = ceil(src/2);
    clone_dist(i,:) = calculate_average_read_distance(T, nMuts, t, src);
end

%%  For each node in each tree store how many reads from each source are
%   below it in tree, and how many reads are in total

clone_homogeny = struct('down', [], 'total', 0);
clone_homogeny = clone_homogeny(ones(N,1));

for i=1:N
    if mod(i,1000)==0, fprintf('%d\n', i); end
    
    [~, W id src Z] =  play_with_clone([tree_dir trees(i).name]);
    T = Z.a.tree;
    t = Z.a.t;
    
    clone_homogeny(i).down = zeros(size(T,1),4);
    clone_homogeny(i).total = T(:,3);
    for s=1:4
        cnt = hist(t(src == s), 1:size(T,1));
        T_ = [T(:,1) cnt(:)];
        T_ = fill_counts(T_);
        clone_homogeny(i).down(:,s) = T_(:,3);
    end
end

%%
clone_src = zeros(N,4);
for i=1:N
    if mod(i,1000)==0, fprintf('%d\n', i); end
    [~, ~, ~, src] =  play_with_clone([tree_dir trees(i).name]);
    clone_src(i,:) = hist(src, 1:4);
end
%% compute an F1 score for each node
% return the max score
max_F1 = zeros(N,4);
for i=1:N
    nNodes = size(clone_homogeny(i).down,1);
    precision = clone_homogeny(i).down ./ clone_homogeny(i).total(:,ones(1,4));
    recall = clone_homogeny(i).down ./ clone_homogeny(i).down(ones(nNodes,1), :);
    %F1 = 2*precision .* recall ./ (precision + recall);
    F1 = 2./ (1./precision + 1./recall);
    max_F1(i,:) = max(F1);
end


%% label the nodes according to their sources
[~,labels] = max(clone_src > 0, [], 2);
ix = sum(clone_src > 0, 2) > 1;
labels(ix) = 0;
cur_label = 5;
for i=1:4
    for j=(i+1):4
        labels(clone_src(:,i)>1 & clone_src(:,j)>1) = cur_label;
        cur_label = cur_label+1;
    end
end
ix = sum(clone_src > 1, 2) > 2;
labels(ix) = cur_label;
max_label = cur_label;
names = {'1', '2', '3', '4', '12', '13', '14', '23', '24', '34', 'triplets'};
%%
clone_dist_ = clone_dist(:,1:4) ./ clone_dist(:,5*ones(1,4));
%% show for every source, a histogram over the pairwise distances
% among the reads of this source in all the shared clones that have this
% source.
close all
figure;
ix = find(sum(clone_src>1, 2) > 1);
x = reshape(clone_dist_(ix,:), [], 1);
labels = repmat(1:4, length(ix),1);
labels = labels(:);
iy = ~isnan(x);
x = x(iy); labels = labels(iy);

%[~, ~, P] = ...
    hist_with_labels(x, labels, 'source-normalized avg read distance / tree avg read distance', 'clones', 0.1, 1, 'separate');            

xlim([-0.2 10]);
ylim([0 40]);

%% scatter plot for the F1 scores of all clones from each pair of sources
figure;
junk = clone_src > 1;
[codes,~,doubles] = get_tuple_codes();
N = zeros(length(doubles),1);
for i=1:length(doubles)    
    src = find(codes(doubles(i),:));
    iy = (junk(:,src(1)) & junk(:,src(2)) );
    N(i) = sum(iy);
    subplot(length(doubles)/2,2,i);
    plot(clone_dist_(iy,src(1)), clone_dist_(iy,src(2)), 'x');
    title(sprintf('sources %d%d', src(1), src(2)));
     xlim([0 5])
     ylim([0 5]);
end
    
%% Pairwise read distances but now not based on tree at all
fasta_file_name = '/afs/cs/u/joni/scratch/data/lymph/t5.fa';
load([fasta_file_name '.mat']);      
data = correct_PCR_artifacts(data);
%%
intra_read = zeros(N,26);
for i=1:N
    if mod(i,1000)==0, fprintf('%d\n', i); end
    [dg dr Dsrc] = analyze_intra_read_similarity([tree_dir trees(i).name], data);
    intra_read(i,1:10) =  dg(:);
    intra_read(i,11:20) = dr(:);
    intra_read(i,21:26) = Dsrc(:);
end



%%
chosen = zeros(N,1);
iy = find(sum(clone_src > 1, 2) > 1);
for i=iy'
    close all hidden
    i
    [clone, G, id, src, Z] = play_with_clone([tree_dir trees(i).name]);    
    clone = clone(2:end,:);

    map('ACGTN') = 1:5;
    priors = struct('is_codon', false, 'is_decay', false, 'pi', map(G.Sequence));
    priors.R = generate_sticky_prior(1e-7, 4, 6300000);
    priors.nClasses = 1;
    priors = init_priors(priors);

    it = 1; ll__ = Z.ll;
    for it=1:size(ll__,2)
        r = Z.cur(it+1);
        if mod(it,1000) == 0, fprintf('%d ', it); end
        a = collapse_edges(convert_phylo_tree_to_mutation_tree(r));
        T = a.tree2;
        [~, ll__(3,it)] = MH_rates(T, r.F, r.rates, 0);    
        [ll__(4:5,it),~] = MH_mutation_parameters(T, a.sequences, a.t, r.mut_model, priors);
    end
    [~,k] = max(sum(ll__));
    chosen(i) = k;
    
    ix = [find(G(1).Header=='_', 6, 'last') length(G(1).Header)+1];
    n = diff(ix(5:7))-1;
    str = sprintf('  %s\nN1/2: [%d %d]; ', trees(i).name, n(1), n(2));
    str = [str sprintf(' eaten: [%s]', G(1).Header(ix(1)+1 : ix(5)-1))];
    str = [str sprintf('\nread2germ: [%.1f %.1f %.1f %.1f %.1f]', intra_read(i,1:5))];
    str = [str sprintf('\nread2read: [%.1f %.1f %.1f %.1f %.1f]', intra_read(i,11:15))];
    str(str == '_') = ' ';
    
    
    a = Z.cur(k+1); %best;
    a = convert_phylo_tree_to_mutation_tree(a);
    a = collapse_edges(a);
    h = visualize_tree(a.tree, a.sequences, [a.t src], clone,4);
    if ~isempty(h)
        h.Nodes(1).Label = [h.Nodes(1).Label sprintf(' iteration %d', k) str];
        jpg_file_name = [tree_dir trees(i).name(1:end-3) '_titled.jpg'];
        saveas(h.hgAxes, jpg_file_name); 
    end        
end
%save([tree_dir sprintf('single_point_%d', i)], 'chosen');



