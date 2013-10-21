function collect_statistics_on_trees()
%%  Loading prerequisites
clear
cd ~/JVL/src/phylo/Fire
addpath('.');
addpath('../../phylo/');
addpath('../../DPtrees/');
addpath('../../VDJ/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));
%tree_dir = '/afs/cs/u/joni/scratch/data/lymph/t5.fa_results_take_13_clones_greater_than_two_reads/';
tree_dir = '/afs/cs/u/joni/scratch/data/lymph/t5.fa_results/clones_not_cleaned/';
trees = dir([tree_dir '/*clone*.fa']);


%% Init parameters
N = length(trees);
clone_ids =   zeros(N,5); % extract the parameters embedded in the filename
clone_stats = zeros(N,10);% For each clone: no. reads and "eaten".
clone_src  =   zeros(N,24); % How many reads from each source/rep
clone_src_ =   zeros(N,4); % How many reads from each source
clone_hash =  zeros(N,8); % a checksum that can detect clones with identical composition of reads
clone_n1dn2 = cell(N,3);  % N1, N2 and D segments of each clone.
M = 3e5;
inds = sparse(M,N, false);  % will which read is in which clone

%hash = @(x) mod(sum(cell2mat(x)-'0'), 100);
hash = @(x) sum(x); % the has function
eaten = zeros(1,6);


%% store all the winning states in one struct.

Z = load([tree_dir trees(1).name(1:end-3) '.mat']);           
states = Z.a;
states(ones(1,N),:);

for i=1:N
    if mod(i,1000)==0, fprintf('%d ', i); end
    W = fastaread([tree_dir trees(i).name]);
    Z = load([tree_dir trees(i).name(1:end-3) '.mat']);           
    if true
        [clone_fasta{i} states(i)] = absorb_leaves(W, Z.a);
    else
        clone_fasta{i}  = W;
        states(i) = Z.a;
    end
end
fprintf('\n');


%% fill in clone_ids, clone_stats, clone_hash, clone_n1dn2, clone_src, inds
% For each tree
err = [];
for i=1:N
    if mod(i,1000)==0, fprintf('%d\n', i); end
    [~, G, read_ix, src_, rep_] =  play_with_clone(clone_fasta{i});
    
    clone_ids(i,:) = sscanf(trees(i).name', 'patient_%d_V_%d_J_%d_len_%d_clone_%d')';
    src = sub2ind([4 6], src_, rep_); % source + replicate information
    clone_src(i,:) =  hist(src, 1:24);
    clone_src_(i,:) = hist(src_,1:4); % just source
    inds(read_ix,i) = true;
    
    clone_stats(i,1) = length(src_);
    clone_stats(i,3:8) = G.eaten;
    clone_stats(i,9:10) = G.lastV_firstJ;
    clone_n1dn2{i,1} = G.n1; 
    if ~isempty(G.n2), clone_n1dn2{i,2} = G.n2; end
    clone_n1dn2{i,3} = G.d; 
end
assert(isempty(err));


%%  throw away clones entirely from one replica - not relevant anymore.
% del = [sum(clone_src,2) == max(clone_src, [],2)];
% del = find(del)';
% fprintf('%d single replica clones are deleted!\n', length(del));
% 
% %%
% I = find(~ismember(1:N, del)); % clones that remain
% clone_src = clone_src(I,:);
% clone_src_ = clone_src_(I,:);
% clone_ids = clone_ids(I,:);
% clone_stats = clone_stats(I,:);
% clone_n1dn2 = clone_n1dn2(I,:);
% trees = trees(I);
% N = length(trees);
% inds = inds(:,I);

%% mark for deletion clones contained in other clones
 
%[~, I, clone_stats(:,2)] = unique(inds', 'rows');
C = double(inds)'*double(inds); % How many reads are shared between each pair of clones 

% extract the diagonal (how many reads are in each clone), and zero it
sum_inds = diag(C);             
C(sub2ind([N N], 1:N,1:N)) = 0;

% Work on pairs that have any overlap
[I1 I2] = find(C > 0); 
ix = sub2ind([N N], I1,I2);

% Each pair is counted twice, so zero one of those times
C(ix(I1<I2)) = 0;
[I1 I2] = find(C > 0);  
ix = sub2ind([N N], I1,I2);

% Calculate matrix D: Jaccard index for each pair of clones
%   intersection/union
union_sets = sum_inds(I1) + sum_inds(I2) - C(ix); 
D = C;
D(ix) = C(ix)./union_sets; 

% If a pair has Jaccard index > threshold, eliminate smaller clone
overlap_pairs = find(D(ix)>0.65);
del = [];
for i=overlap_pairs'
    fprintf('%s (%d) and\n%s (%d) --> %.1f%% overlap.', trees(I1(i)).name, ...
        sum_inds(I1(i)), trees(I2(i)).name, sum_inds(I2(i)), 100*D(I1(i), I2(i))); 
    if sum_inds(I1(i)) <= sum_inds(I2(i))
        del = [del I1(i)];
        fprintf('  --> Erasing first.\n\n');
    else 
        del = [del I2(i)];
        fprintf('  --> Erasing second.\n\n');
    end
end
del = unique(del);             % clones to delete
fprintf('%d pairs of clones with overlap.\n', length(del));


%% Mark for deletion clones with reads of patient 1's spleen
%locate clones
ix = find(clone_ids(:,1) == 1);
iy = find(clone_src_(ix,2) > 0);
del = unique([del ix(iy)']);



%%  move deleted clones files - no need!
% Build script to move the files corresponding to deleted clones to
%  "redundant" directory
% fid = fopen([tree_dir 'move_redundant.sh'], 'w');
% fprintf(fid, ['mkdir ' tree_dir '/redundant\n']);
% for i=del
%     str = trees(i).name(1:end-3);
%     fprintf(fid, 'mv %s%s.* %sredundant\n', tree_dir, str, tree_dir);
% end
% fclose(fid);
%% Reinitialize variables with remaining clones

fprintf('Deleting %d clones.\n', length(del));
I = find(~ismember(1:N, del)); % clones that remain
clone_src = clone_src(I,:);
clone_src_ = clone_src_(I,:);
clone_ids = clone_ids(I,:);
clone_stats = clone_stats(I,:);
clone_n1dn2 = clone_n1dn2(I,:);
trees = trees(I);
N = length(trees);
inds = inds(:,I);
states = states(I); 
clone_fasta = clone_fasta(I);



%%  Load the tree for each clone, fill in clone_mat.
% distance from root
clone_mat = struct('height', 0, 'nNodes', 0, 'medMut', 0, 'dist_to_root', 0, 'maxMut', 0, 'len', 0, 'med_depth', 0, 'med_depth_leaves', 0,...
    'med_depth_nonempty', 0, 'avg_node_distance', 0);
clone_mat = clone_mat(ones(1,N),:);

for i=1:N      
    if mod(i,1000)==0, fprintf('%d\n', i); end
    st = states(i);

    assert(st.tree(1,1) == 0);
    assert(st.tree(2,1) == 1);
    [nMuts mut] = annotate_mutations_on_tree(st.tree(:,1), st.sequences);        
    % recall the function annotate_mutations does not count mismatches with
    % gaps as a mutation.  

    % distance to germline
    clone_mat(i).dist_to_root = nMuts(2); % exclude junction mutations
    nMuts(2) = 0;

%    clone_dist(i,:) = calculate_average_read_distance(st.tree, nMuts, st.t, src);
    clone_mat(i).avg_node_distance = calculate_average_read_distance(st.tree, nMuts, []);

    mut = sum(mut ~= -2, 1);
    
    % depth of tree (not including germline to root of clone)
    height = fill_up([st.tree(:,1) nMuts], @max, false);
    depth  = fill_down([st.tree(:,1) nMuts], @sum);

    clone_mat(i).height = height(1,3);
    clone_mat(i).avg_height = mean(height(2:end,3));
    clone_mat(i).med_depth = mean(depth(2:end,3));
    clone_mat(i).med_depth_leaves = mean(depth(st.tree(:,2)==st.tree(:,3),3));
    clone_mat(i).med_depth_nonempty = mean(depth(st.tree(:,2)>0,3));
    assert(max(depth(:,3)) == height(1,3));
    
    clone_mat(i).leaves_frac = sum(st.tree(:,2)==st.tree(:,3)) / (size(st.tree,1)-1);

    len = size(st.sequences,2);
    clone_mat(i).len = len;
    %  if reading-frame corrent:  length + whatever was trimmed mod 3 = 1
    clone_mat(i).RF = mod(len + clone_stats(i,3) + clone_stats(i,8) - 1, 3);
    clone_mat(i).nNodes = size(st.sequences,1);     
    seg = [clone_stats(i,9) clone_stats(i,10)];
    clone_mat(i).mut =  [sum(mut(1:seg(1))) ...
                         sum(mut(seg(1)+1:seg(2)-1)) ...
                         sum(mut(seg(2):end)) ];
    clone_mat(i).mut =  clone_mat(i).mut ./ [seg(1) seg(2)-seg(1)-1 len-seg(2)+1];
    
                     
    % median number of mutations
    nMuts(nMuts == 0) = [];
    if isempty(nMuts), nMuts = 0; end
    clone_mat(i).medMut = median(nMuts);    
    clone_mat(i).maxMut = max(nMuts);   


end
    % length of junctions  

len_n1 = cellfun(@length, clone_n1dn2(:,1));
len_n2 = cellfun(@length, clone_n1dn2(:,2));
for i=1:N, 
    clone_mat(i).len_n1 = len_n1(i);
    clone_mat(i).len_n2 = len_n2(i);
end


%% More clone stats
clone_intra_read = zeros(N,20);
for i=1:N
    if mod(i,1000)==0, fprintf('%d\n', i); end
    [dg dr] = analyze_intra_read_similarity(clone_fasta{i}, true, false);
    clone_intra_read(i,1:10) =  dg(:);
    clone_intra_read(i,11:20) = dr(:);
%    clone_intra_read(i,21:26) = Dsrc(:);
end
 
  
%% post process mutation model take counts from all trees
I = (eye(4) ~= 1);

rep = load_repertoire('ihmmune_collapsed');
A = get_IGH_alignments(rep, 'V');

clone_mutations = struct('silent', -1,  'third_base', -1, 'single_base', -1, 'VDJ', [], 'VDJ_per_base', [], ...
    'hotspots', [], 'hotspotsV', [], 'root', [], 'silentspots', [], 'nonsilentspots', [], 'germ_silentspots', [], ...
    'germ_nonsilentspots', [], ...
    'al_hotspots', zeros(1,size(A,2), 'int8'), 'al_germ_hotspots', zeros(1,size(A,2), 'int8'));
clone_mutations = clone_mutations(ones(1,N),1);

counts = zeros(4,4); 
problem = [];
for i=1:N
    if mod(i,1000)==0, fprintf('%d ', i); end  
    v = clone_ids(i,2);
    
    T = states(i).tree(:,1);
    sequences_c = states(i).sequences(  3:end ,:);
    sequences_p = states(i).sequences(T(3:end),:);
    x = (sequences_p-1)*4 + sequences_c;
    y = sum(I(x),1); % y(l) is the #mutations in location l    
    seg = clone_stats(i,9:10);
    clone_mutations(i).VDJ = [sum(y(1:seg(1))) sum(y(seg(1)+1:seg(2)-1)) sum(y(seg(2):end))];    
    clone_mutations(i).VDJ_per_base = clone_mutations(i).VDJ ./ ...
                                      [seg(1) seg(2)-seg(1)-1 length(y)-seg(2)+1];        
    counts = counts + reshape(histc(x(:), 1:numel(counts)), size(counts));
    
    % show hotspots normalized by lengh of whole sequence
    len = length(y);
    bins = ceil(100* (1:len)/ len );
    J = false(len,100);
    J(sub2ind(size(J), 1:len, bins)) = true;
    clone_mutations(i).hotspots= int8(y*J);        

    
    % show hotspots on V using IMGT alignments    %%%%%
    y = [zeros(1,clone_stats(i,3)) y(1:clone_stats(i,9)) zeros(1,clone_stats(i,4))];
    try
        assert(length(y) == length(rep.V(v).Sequence));
        assert(length(y) == sum(A(v,:)));
        clone_mutations(i).al_hotspots(A(v,:)) = y;
    catch 
        problem = [problem i];
    end       
    %%%%%%%%%%%%%%%%%%%%%%                        %%%%%
    
    % show hotspots  normalized by length in root of clone compared to
    % germline
    y = double((states(i).sequences(1, :) ~= states(i).sequences(2, :)) & (states(i).sequences(1, :) ~= 5));
    clone_mutations(i).root = int8(y*J);                
    
    % show hotspots normalized by lengh of V region
    len = seg(1);
    bins = ceil(50* (1:len)/ len );
    J = false(len,50);
    J(sub2ind(size(J), 1:len, bins)) = true;
    clone_mutations(i).hotspotsV = int8(y(1:seg)*J); 

    % show hotspots on V using IMGT alignments    %%%%%
    y = [zeros(1,clone_stats(i,3)) y(1:clone_stats(i,9)) zeros(1,clone_stats(i,4))];
    clone_mutations(i).al_germ_hotspots(A(v,:)) = y;
    
end
fprintf('\n');

% show transition matrix obtained from counts
R = generate_stochastic_matrix(counts, 1, false);


%% post process for silent/non_silent

% find all the clones that are in frame
% junk = [ mod(clone_stats(:,3),3) [clone_mat.RF]' ];
% ix = find( junk(:,1) == 0 & junk(:,2) == 0)';
[~, I, J, K] = silent_map();

counts = zeros(64,64); 
for i= find([clone_mat.RF] == 0)
    if mod(i,1000)==0, fprintf('%d ', i); end 
    v = clone_ids(i,2);

    T = states(i).tree(:,1);
    
    % trim until first codon base
    Vtrim = mod(clone_stats(i,3),3);
    sequences   = states(i).sequences(:, Vtrim+1:end);

    % convert to codons and trim from the right
    sequences  = seqs2codons(sequences(:,1:end-mod(size(sequences,2),3)));
    
    % fill in junctions of germline with root of tree 
    sequences(1,sequences(1,:) == 65) = sequences(2,sequences(1,:) == 65);
    assert(max(sequences(:)) < 65);
    sequences_p = sequences(T(3:end),:);
    sequences_c = sequences(  3:end, :);
    x = (sequences_p-1)*64 + sequences_c;
    clone_mutations(i).single_base = sum(J(x(:))); % single-base mutation
    clone_mutations(i).silent =      sum(I(x(:))); % single-base silent mutation
    clone_mutations(i).third_base =  sum(K(x(:))); % single-base third base mutation    
    counts = counts + reshape(histc(x(:), 1:numel(counts)), size(counts));
    
    % show silent/nonsilent mutations normalized by lengh of whole sequence
    len = size(sequences,2);
    bins = ceil(30* (1:len)/ len );
    F = false(len,30);
    F(sub2ind(size(F), 1:len, bins)) = true;
    
    y = sum(I(x),1);    % silent
    clone_mutations(i).silentspots = int8(y*F);        
    y = sum(J(x)-I(x),1);    % non-silent
    clone_mutations(i).nonsilentspots= int8(y*F);        

    x = (sequences(1,:)-1)*64 + sequences(2,:);
    y = sum(I(x),1);    % silent
    clone_mutations(i).germ_silentspots = int8(y*F);        
    y = sum(J(x)-I(x),1);  % non-silent
    clone_mutations(i).germ_nonsilentspots= int8(y*F);   
    
end
fprintf('\n');

% print fraction of silent mutations
counts_ = counts;
counts_(~J) = 0;
fprintf('fraction of silent mutations:\n');
sum(counts_(I))/ sum(counts_(J))
fprintf('fraction of mutations in 3 base:\n');
sum(counts_(K))/ sum(counts_(J))


%% post process for silent/non_silent JUST IN V REGION

rep = load_repertoire('ihmmune_collapsed');
[~, I, J, K] = silent_map();
A = get_IGH_alignments(rep, 'V');

for i=1:N
    if mod(i,1000)==0, fprintf('%d ', i); end 
    v = clone_ids(i,2);

    T = states(i).tree(:,1);

    nNodes = size(states(i).sequences,1);

    % the parts of the sequences that were not trimmed/eaten - bases that
    % can have mutations
    support    = [zeros(1,clone_stats(i,3))  ones(1,clone_stats(i,9)) zeros(1, clone_stats(i,4))];
    al_support = zeros(1,size(A,2));
    al_support(1, A(v,:)) = support; % aligning to reference V
    al_support = al_support(1:end-mod(size(al_support,2),3));
    al_support = min(reshape(al_support, 3, []));
    clone_mutations(i).V.support = al_support;

    
        
    sequences  = [5*ones(nNodes,clone_stats(i,3))  states(i).sequences(:, 1:clone_stats(i,9)) 5*ones(nNodes, clone_stats(i,4))];
    assert(size(sequences,2) == length(rep.V(v).Sequence));
    Al_V = 5*ones(nNodes,size(A,2));
    Al_V(:, A(v,:)) = sequences;
    sequences  = seqs2codons(Al_V(:,1:end-mod(size(Al_V,2),3)));
    assert(all(unique(sum(sequences == 65,1)) == [0 nNodes]));
    sequences(sequences == 65) = 1;
    sequences_p = sequences(T(3:end),:);
    sequences_c = sequences(  3:end, :);
    
    % TODO: Merge the following two blocks
    x = (sequences_p-1)*64 + sequences_c;               
    y = sum(I(x),1);    % silent
    clone_mutations(i).V.silentspots = int8(y);        
    y = sum(J(x)-I(x),1);    % non-silent
    clone_mutations(i).V.nonsilentspots= int8(y);        

    x = (sequences(1,:)-1)*64 + sequences(2,:);
    y = sum(I(x),1);    % silent
    clone_mutations(i).V.germ_silentspots = int8(y);        
    y = sum(J(x)-I(x),1);  % non-silent
    clone_mutations(i).V.germ_nonsilentspots= int8(y);   
    
end
fprintf('\n');


%% post process for silent/non_silent JUST IN J REGION
% TODO: Merge this code with the V code.  

rep = load_repertoire('ihmmune_collapsed');
[~, I, J, K] = silent_map();
A = get_IGH_alignments(rep, 'J');

for i=1:N
    if mod(i,1000)==0, fprintf('%d ', i); end 
    j = clone_ids(i,3);

    T = states(i).tree(:,1);

    support    = [zeros(1,clone_stats(i,7))  ones(1,size(states(i).sequences,2)-clone_stats(i,10)+1) zeros(1, clone_stats(i,8))];
    J_support = zeros(1,size(A,2));
    J_support(1, A(j,:)) = support;
    J_support = J_support(1:end-mod(size(J_support,2),3));
    J_support = min(reshape(J_support, 3, []));    
    clone_mutations(i).J.support = J_support;
    

    nNodes = size(states(i).sequences,1);
    sequences  = [5*ones(nNodes,clone_stats(i,7))  states(i).sequences(:, clone_stats(i,10):end) 5*ones(nNodes, clone_stats(i,8))];
    assert(size(sequences,2) == length(rep.J(j).Sequence));
    Al_J = 5*ones(nNodes,size(A,2));
    Al_J(:, A(j,:)) = sequences;
    sequences  = seqs2codons(Al_J(:,1:end-mod(size(Al_J,2),3)));
    assert(all(unique(sum(sequences == 65,1)) == [0 nNodes]));
    sequences(sequences == 65) = 1;
    sequences_p = sequences(T(3:end),:);
    sequences_c = sequences(  3:end, :);

    % TODO: Merge the following two blocks
    x = (sequences_p-1)*64 + sequences_c;               
    y = sum(I(x),1);    % silent
    clone_mutations(i).J.silentspots = int8(y);        
    y = sum(J(x)-I(x),1);    % non-silent
    clone_mutations(i).J.nonsilentspots= int8(y);        

    x = (sequences(1,:)-1)*64 + sequences(2,:);
    y = sum(I(x),1);    % silent
    clone_mutations(i).J.germ_silentspots = int8(y);        
    y = sum(J(x)-I(x),1);  % non-silent
    clone_mutations(i).J.germ_nonsilentspots= int8(y);   
    

end
fprintf('\n');

%%  Bushy-Chainy and other triplet based analysis
clone_triplets =  struct('chainy', NaN, 'bushy', NaN);
clone_triplets = clone_triplets(ones(1,N), 1);

for i=1:N
    if mod(i,1000)==0, fprintf('%d ', i); end
    % prepare list of triplets
    if clone_stats(i,1) < 3,  continue; end
    
    st = states(i);
    clone =  play_with_clone(clone_fasta{i});
            
    if size(clone,1) == 3
        triplets = [1 2 3];
    elseif size(clone,1) == 4 
        triplets = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
    elseif size(clone,1) == 5
        triplets = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
        triplets = [triplets ; 1+triplets(2:end,:) ; 1 4 5; 1 2 5; 1 3 5];
    else
        triplets = ceil(size(clone,1)*rand(100,1));
        triplets = [triplets ceil((size(clone,1)-1)*rand(100,1))];
        triplets = [triplets ceil((size(clone,1)-2)*rand(100,1))];
        ix = triplets(:,2) >= triplets(:,1);
        triplets(ix,2) = triplets(ix,2)+1;
        ix = triplets(:,3) >= triplets(:,2);
        triplets(ix,3) = triplets(ix,3)+1;        
    end
    
    % calcualte D, distance matrix between reads. Exclude mismatches with gaps.
    B = pdist(double(clone), 'hamming')*size(clone,2);
    R = pdist(double(clone == 5), 'hamming')*size(clone,2);
    D = B-0.5*R;    % the 0.5 is a necessary hack to make sure that the 
                    % triangle inequility still holds
    D = squareform(D);
    
    % given 3 reads, x, y, z
    % For each read, find if it is in the path between the other two by
    % computing for each read y U(y) = D(x,y)+D(y,z)-D(x,z).    
    V = [1 1 -1; 1 -1 1; -1 1 1];
    ix =     sub2ind(size(D), triplets(:,1), triplets(:,2));
    ix = [ix sub2ind(size(D), triplets(:,2), triplets(:,3))];
    ix = [ix sub2ind(size(D), triplets(:,1), triplets(:,3))];
    U = D(ix)*V;    
    
    % the bushiness of a triplet x y z is min(U(x) U(y) U(z))
    % the bushiness of a clone is the average bushiness of triplets
    clone_triplets(i).bushy = 0.5*mean(min(U,[],2));
    
    % another approach:  Find how often one read is between the other two
    % on the tree.
    
    % Construct D.  D(i,j) = true <==> j is ancestor of i,
    T = st.tree;
    M = size(T,1);
    [anc, o, d] = order_tree(st.tree);
    D = false(M,M);
    for j=1:length(o)
        D(o(j), anc(j,1:d(j))) = true;
    end
    
    
    cnt = 0;
    for j=1:size(triplets,1)
        % trp indicates the tree nodes of the reads in the current triplet
        trp = false(1,M);    
        trp(st.t(triplets(j,:))) = true;
        
        % if all the nodes in trp appear as ancestors anywhere ==> cnt++
        ix = sum(trp(ones(M,1),:) <= D,2);
        cnt = cnt + any(ix == M);
    end
%    cnt/size(triplets,1)
    
    clone_triplets(i).chainy = cnt/size(triplets,1);
    
end
fprintf('\n');



%% save state
save([tree_dir 'all_clones.mat'], 'N', 'clone_ids', ...
    'clone_mat', 'clone_n1dn2', 'clone_src', 'clone_src_', 'clone_stats', ...
    'inds', 'trees', 'clone_intra_read', 'clone_mutations', 'clone_triplets');

%%
save([tree_dir 'states.mat'], 'states', 'clone_fasta');



