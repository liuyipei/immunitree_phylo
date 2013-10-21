function b = load_newick_tree(filename, leaves_with_one_read)

if nargin<2, leaves_with_one_read = false; end
h = phytreeread(filename);
[MATRIX, ID] = getmatrix(h);

N = (size(MATRIX,1)+1)/2;
ix = zeros(N,2); 
for i=1:N
    ix(i,:) = sscanf(ID{i}, 'id_%d_node_id_%d'); 
end
[~, order] = sort(ix(:,1));

[~, parents] = max(MATRIX);
b = [];
b.tree = parents';
assert(b.tree(end) == 1);
b.tree(1:N) = b.tree(order);
if leaves_with_one_read
    b.t = 1:N;
else
    b.tree = b.tree-N;
    b.t = b.tree(1:N);
    b.tree = b.tree(N+1:end);
end

M = size(b.tree,1);
b.tree(end) = 0;
b.tree = build_tree_counts(b.tree, b.t);
b.labels = ix(order,2);
b.sequences = ones(M,1);
b.tree2 = [b.tree(:,1) zeros(M,1), 10*ones(M,1)];

end


function T = build_tree_counts(T,t)
    h = hist(t, size(T,1))';
    T = [T h];
    T = fill_counts(T);
end

function test()
%%
X = fastaread('/afs/cs/u/joni/JVL/src/phylo/tandy.fa');
map('ACGTN') = 1:5;
reads = int16(map(cell2mat({X.Sequence}')));
b = load_newick_tree('/afs/cs/u/joni/scratch/software/FastTree_2.0/TEST_DEC_29_2011.newick')
% [names{1001:end}] = deal([]);
b.sequences = ones(size(b.tree,1), size(reads,2));

mut_model.NT = generate_sticky_prior(1e-3, 4, 1);
mut_model.R = [generate_sticky_prior(1e-3, 4, 1); ones(1,4)]; 

global greedy
greedy = -1;
b.sequences = gibbs_sequences_of_nodes_with_evidence(b.tree2,b.sequences, 10, ...
    mut_model.NT(:,:,ones(1,size(reads,2))), [], b.t, mut_model.R, reads);

a = convert_phylo_tree_to_mutation_tree(b);
a = collapse_edges(a);
a = order_nodes(a)

view_tree(a.tree, a.sequences, 1, [a.t b.labels]);


end