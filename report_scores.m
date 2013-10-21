function [tendency, rand_index] = report_scores(s, a, labels)

a = convert_phylo_tree_to_mutation_tree(a);
D_ = node2node_dist(a.tree, a.sequences);
B_ = squareform(D_(a.t,a.t));


if 0 %%~isstruct(s)
    reads = s;
    s = [];
    [s.tree s.sequences s.t] = initialize_chain(reads);
end

if isstruct(s)
    s = convert_phylo_tree_to_mutation_tree(s);
    D  = node2node_dist(s.tree, s.sequences);
    B  = squareform(D(s.t,s.t));    
    if ~exist('labels', 'var'), labels = s.t; end
else % s is all the reads
    [sequences, ~, t] = unique(s, 'rows'); 
    D = squareform(pdist(double(sequences), 'hamming')*size(sequences,2));
    B  = squareform(D(t,t));
end

% Compute RAND index
if exist('labels', 'var')
    rand_index = confusion_matrix(labels, a.t)
end

% consider the average distance between two reads.
max_dist = max(max(D(:)), max(max(D_(:))));
Z = B_-B;

N = length(a.t);
total_pairs = N*(N-1)/2;

were_together = sum(B==0);

% how many pairs of reads where in distance 0 before?  Where are they now?
still_together = sum(B == 0 & B_ == 0);

% no. of times dist hasn't changed.
same_dist = sum(Z == 0);


avg_dist_change = sum(abs(Z))/(total_pairs - same_dist);
tendency = sum(Z)/(total_pairs - same_dist);
if tendency > 0, trend_str = 'longer'; else trend_str = 'shorter'; end


str = sprintf('There is a total of %d total pairs (no of reads = %d).\n', total_pairs, N);
str = [str sprintf('Of these pairs, %.1f%% are in the same class in truth.\n', ...
    100*were_together/total_pairs)];
str = [str sprintf('Of these pairs, %.1f%% are still together now, which is %.1f%% of the total pairs.\n',...
    100*still_together/were_together, 100*still_together/total_pairs)];
str = [str sprintf('%.1f%% of all pairs have same tree distance before and after.\n',...
    100*same_dist/total_pairs)];
str = [str sprintf('Those pairs that changed their distance, changed it by %.1f edges on average.\n',...
    avg_dist_change)];
str = [str sprintf('Overall, the distances are now *%s* by %.1f.\n', trend_str, abs(tendency))];

if nargout == 0
    hist(Z, -max_dist:max_dist);
    title(str);
%    title(sprintf('change in read-pairs distance\nOn average, the distances are now *%s* by %.1f.', trend_str, abs(tendency)));
    xlabel('distance change');
    ylabel('no. of pairs');
end
%stats(k,j,:) = [total_pairs were_together still_together same_dist avg_dist_change tendency rand_index];

end



function D = node2node_dist(T, seqs)
    nMuts = annotate_mutations_on_tree(T, seqs);
    assert( (nMuts(1) == 0) && (sum(nMuts==0) == 1) );
    DG = get_graph_from_tree(T, nMuts, true);
    D = graphallshortestpaths(DG, 'directed', false);
end

