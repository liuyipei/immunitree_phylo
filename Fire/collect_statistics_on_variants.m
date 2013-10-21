function collect_statistics_on_variants()
%%  Loading prerequisites
clear
cd ~/JVL/src/phylo/Fire
addpath('.');
addpath('../../phylo/');
addpath('../../DPtrees/');
addpath('../../VDJ/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));
tree_dir = '/afs/cs/u/joni/scratch/data/lymph/t5.fa_results/nonclones/';
%trees = dir([tree_dir 'patient*.fa']);
trees = dir([tree_dir '*.fa']);

% for a new set run sep_file.sh

%% Init parameters
N = 200000; %length(trees);
clone_ids =   zeros(N,3); % extract the parameters embedded in the filename
clone_stats = zeros(N,11);% For each clone: no. reads and "eaten".
clone_src  =   zeros(N,24); % How many reads from each source/rep
clone_src_ =   zeros(N,4); % How many reads from each source
clone_n1dn2 = cell(N,3);  % N1, N2 and D segments of each clone.
clone_fasta = cell(N,1);  % N1, N2 and D segments of each clone.
M = 3e5;
inds = sparse(M,N, false);  % will which read is in which clone

%hash = @(x) mod(sum(cell2mat(x)-'0'), 100);
hash = @(x) sum(x); % the has function
eaten = zeros(1,6);
%% fill in clone_ids, clone_stats, clone_n1dn2, clone_src, inds - erase me
% For each tree
clone_src   = zeros(N,240); % How many reads from each source/rep
clone_src_  = zeros(N,40); % How many reads from each source
clone_src__ = zeros(N,10);
i = 0;
for k = 1:length(trees) 
    if mod(k,100)==0, fprintf('%d\n', i); end
    patient_v_j = sscanf(trees(k).name', 'V_%d_J_%d')';
    W = fastaread([tree_dir trees(k).name]);
    J = strmatch('germline', {W.Header});
    J = [J; length(W)+1];
    for j=1:length(J)-1        
        i = i+1;
        clone_fasta{i} = W(J(j):(J(j+1)-1));
        [~, G, read_ix, src_, rep_] =  play_with_clone(clone_fasta{i});

        clone_ids(i,2:3) = patient_v_j;
        src = sub2ind([40 6], src_, rep_); % source + replicate information
        clone_src(i,:) =  hist(src, 1:240);
        clone_src_(i,:) = hist(src_,1:40); % just source
        inds(read_ix,i) = true;

        clone_stats(i,1) = length(src_);
        clone_stats(i,2) = length(G.Sequence);
        clone_stats(i,3:8) = G.eaten;
        clone_stats(i,9:10) = G.lastV_firstJ;
        clone_stats(i,11) = mod(clone_stats(i,2) + clone_stats(i,3) + clone_stats(i,8) - 1, 3);
        clone_n1dn2{i,1} = G.n1; 
        if ~isempty(G.n2), clone_n1dn2{i,2} = G.n2; end
        clone_n1dn2{i,3} = G.d; 
    end
end
clone_src_(:,2) = 0;
clone_src(:,2:40:end) = 0;
N = i;

%% fill in clone_ids, clone_stats, clone_n1dn2, clone_src, inds
% For each tree
i = 0;
for k=1:length(trees) 
    if mod(k,100)==0, fprintf('%d\n', i); end
    patient_v_j= sscanf(trees(k).name', 'patient_%d_V_%d_J_%d')';
    W = fastaread([tree_dir trees(k).name]);
    J = strmatch('germline', {W.Header});
    J = [J; length(W)+1];
    for j=1:length(J)-1        
        i = i+1;
        clone_fasta{i} = W(J(j):(J(j+1)-1));
        [~, G, read_ix, src_, rep_] =  play_with_clone(clone_fasta{i});

        clone_ids(i,:) = patient_v_j;
        src = sub2ind([4 6], src_, rep_); % source + replicate information
        clone_src(i,:) =  hist(src, 1:24);
        clone_src_(i,:) = hist(src_,1:4); % just source
        inds(read_ix,i) = true;

        clone_stats(i,1) = length(src_);
        clone_stats(i,2) = length(G.Sequence);
        clone_stats(i,3:8) = G.eaten;
        clone_stats(i,9:10) = G.lastV_firstJ;
        clone_stats(i,11) = mod(clone_stats(i,2) + clone_stats(i,3) + clone_stats(i,8) - 1, 3);
        clone_n1dn2{i,1} = G.n1; 
        if ~isempty(G.n2), clone_n1dn2{i,2} = G.n2; end
        clone_n1dn2{i,3} = G.d; 
    end
end
N = i;

%%
I = 1:N; 
clone_src = clone_src(I,:);
clone_src_ = clone_src_(I,:);
clone_ids = clone_ids(I,:);
clone_stats = clone_stats(I,:);
clone_n1dn2 = clone_n1dn2(I,:);
clone_fasta = clone_fasta(I,:);
inds = inds(:,I);

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
    str1 = sprintf('patient_%d_V_%d_J_%d',clone_ids(I1(i),:));
    str2 = sprintf('patient_%d_V_%d_J_%d',clone_ids(I2(i),:));    
    fprintf('%s (%d) and\n%s (%d) --> %.1f%% overlap.', str1, ...
        sum_inds(I1(i)), str2, sum_inds(I2(i)), 100*D(I1(i), I2(i))); 
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

%% Reinitialize variables with remaining clones

fprintf('Deleting %d clones.\n', length(del));
I = find(~ismember(1:N, del)); % clones that remain
clone_src = clone_src(I,:);
clone_src_ = clone_src_(I,:);
clone_ids = clone_ids(I,:);
clone_stats = clone_stats(I,:);
clone_n1dn2 = clone_n1dn2(I,:);
clone_fasta = clone_fasta(I,:);
N = length(I);
inds = inds(:,I);


%% More clone stats
clone_intra_read = zeros(N,20);
for i=1:N
    if mod(i,1000)==0, fprintf('%d\n', i); end
    [dg dr] = analyze_intra_read_similarity(clone_fasta{i}, false, true);
    clone_intra_read(i,1:10) =  dg(:);
    clone_intra_read(i,11:20) = dr(:);
%    clone_intra_read(i,21:26) = Dsrc(:);
end


    
    

%% post process mutation model take counts from all trees
I = (eye(4) ~= 1);
rep = load_repertoire('ihmmune_collapsed');
A = get_IGH_alignments(rep, 'V');
B = get_IGH_alignments(rep, 'J');

%load([tree_dir 'fastas.mat'])

structV = struct('germ_hotspots', zeros(1,size(A,2), 'int8'), ...
    'germ_silentspots', [],  'germ_nonsilentspots', []);
structJ = struct('germ_hotspots', zeros(1,size(B,2), 'int8'), ...
    'germ_silentspots', [],  'germ_nonsilentspots', []);

clone_mutations = struct('silent', -1,  'third_base', -1, 'single_base', -1, ...
    'germ_silentspots', [],  'germ_nonsilentspots', [], 'V', structV, 'J', structJ);
clone_mutations = clone_mutations(ones(1,N),1);
%% nucleotide level analysis of both V and J regions

map('ACGTN') = 1:5;
for i=1:N
    if mod(i,1000)==0, fprintf('%d ', i); end  
    v = clone_ids(i,2);
    j = clone_ids(i,3);

    
%    [clone germline id src replicas] =  play_with_clone(clone_fasta{i});
    
    % show hotspots  normalized by length in root of clone compared to
    % germline
    sequences = map(cell2mat({clone_fasta{i}.Sequence}'));
    consensus = mode(sequences(2:end,:), 1);
    root = sequences(1,:);
    y = double( (root ~= consensus) & (root ~= 5));
    z = mysub2ind(4, root, consensus, 1);  % marks all mutations
    z(root == 5) = 0;
    
    % show hotspots on V using IMGT alignments    %%%%%
%     x = [zeros(1,clone_stats(i,3)) y(1:clone_stats(i,9)) zeros(1,clone_stats(i,4))];
%     clone_mutations(i).al_germ_hotspots(A(v,:)) = x;

    % mark transitions on V using IMGT alignments
    x = [zeros(1,clone_stats(i,3)) z(1:clone_stats(i,9)) zeros(1,clone_stats(i,4))];    
    clone_mutations(i).V.germ_hotspots(A(v,:)) = x;

    % mark transitions on J using IMGT alignments
    x = [zeros(1,clone_stats(i,7)) z(clone_stats(i,10):end) zeros(1,clone_stats(i,8))];    
    clone_mutations(i).J.germ_hotspots(B(j,:)) = x;
    
        
if 0 % no one cares about that anymore.
   
    % show hotspots normalized by lengh of whole sequence
    len = length(y);
    bins = ceil(100* (1:len)/ len );
    J = false(len,100);
    J(sub2ind(size(J), 1:len, bins)) = true;
    clone_mutations(i).root = int8(y*J);                
    
    % show hotspots normalized by lengh of V region
    seg = clone_stats(i,9:10);
    len = seg(1);
    bins = ceil(50* (1:len)/ len );
    J = false(len,50);
    J(sub2ind(size(J), 1:len, bins)) = true;
    clone_mutations(i).hotspotsV = int8(y(1:seg)*J); 
end
    
end
fprintf('\n');

%% post process for silent/non_silent JUST IN V REGION
rep = load_repertoire('ihmmune_collapsed');
A = get_IGH_alignments(rep, 'V');

map('ACGTN') = 1:5;

% find all the clones that are in frame
% junk = [ mod(clone_stats(:,3),3) [clone_mat.RF]' ];
% ix = find( junk(:,1) == 0 & junk(:,2) == 0)';
[~, I, J, K] = silent_map();

for i= 1:N
    if mod(i,1000)==0, fprintf('%d ', i); end 
    v = clone_ids(i,2);

    % the parts of the sequences that were not trimmed/eaten - bases that
    % can have mutations
    support    = [zeros(1,clone_stats(i,3))  ones(1,clone_stats(i,9)) zeros(1, clone_stats(i,4))];
    V_support = zeros(1,size(A,2));
    V_support(1, A(v,:)) = support; % aligning to reference V
    V_support = V_support(1:end-mod(size(V_support,2),3));
    V_support = min(reshape(V_support, 3, []));
    clone_mutations(i).V.support = int8(V_support);

    sequences = map(cell2mat({clone_fasta{i}.Sequence}'));
    consensus = mode(sequences(2:end,:),1);
    root = sequences(1,:);
    sequences = [root; consensus];
    
    nNodes = size(sequences,1);
    sequences  = [5*ones(nNodes,clone_stats(i,3))  sequences(:, 1:clone_stats(i,9)) 5*ones(nNodes, clone_stats(i,4))];
    assert(size(sequences,2) == length(rep.V(v).Sequence));
    Al_V = 5*ones(nNodes,size(A,2));
    Al_V(:, A(v,:)) = sequences;
    sequences  = seqs2codons(Al_V(:,1:end-mod(size(Al_V,2),3)));
    sequences(2,sequences(2,:) == 65) = sequences(1,sequences(2,:) == 65);
    assert(all(unique(sum(sequences == 65,1)) == [0 nNodes]));
    sequences(sequences == 65) = 1;

    x = (sequences(1,:)-1)*64 + sequences(2,:);
    y = sum(I(x),1);    % silent
    clone_mutations(i).V.germ_silentspots = int8(y);        
    y = sum(J(x)-I(x),1);  % non-silent
    clone_mutations(i).V.germ_nonsilentspots= int8(y);   
    

end
fprintf('\n');

%% post process for silent/non_silent JUST IN J REGION
rep = load_repertoire('ihmmune_collapsed');
A = get_IGH_alignments(rep, 'J');
map('ACGTN') = 1:5;

% find all the clones that are in frame
[~, I, J, K] = silent_map();

for i= 1:N
    if mod(i,1000)==0, fprintf('%d ', i); end 
    j = clone_ids(i,3);

    support    = [zeros(1,clone_stats(i,7))  ones(1,size(clone_fasta{i}(1).Sequence,2)-clone_stats(i,10)+1) zeros(1, clone_stats(i,8))];
    J_support = zeros(1,size(A,2));
    J_support(1, A(j,:)) = support;
    J_support = J_support(1:end-mod(size(J_support,2),3));
    J_support = min(reshape(J_support, 3, []));    
    clone_mutations(i).J.support = int8(J_support);
    

    sequences = map(cell2mat({clone_fasta{i}.Sequence}'));
    consensus = mode(sequences(2:end,:),1);
    root = sequences(1,:);
    sequences = [root; consensus];
    
    nNodes = size(sequences,1);
    sequences  = [5*ones(nNodes,clone_stats(i,7))  sequences(:, clone_stats(i,10):end) 5*ones(nNodes, clone_stats(i,8))];
    assert(size(sequences,2) == length(rep.J(j).Sequence));
    Al_J = 5*ones(nNodes, size(A,2) );
    Al_J(:, A(j,:)) = sequences;
    sequences  = seqs2codons(Al_J(:,1:end-mod(size(Al_J,2),3)));
    sequences(2,sequences(2,:) == 65) = sequences(1,sequences(2,:) == 65);
    assert(all(unique(sum(sequences == 65,1)) == [0 nNodes]));
    sequences(sequences == 65) = 1;

    x = (sequences(1,:)-1)*64 + sequences(2,:);
    y = sum(I(x),1);    % silent
    clone_mutations(i).J.germ_silentspots = int8(y);        
    y = sum(J(x)-I(x),1);  % non-silent
    clone_mutations(i).J.germ_nonsilentspots= int8(y);      
end
fprintf('\n');

%% save state
save([tree_dir 'all_clones.mat'], 'N', 'clone_ids', ...
    'clone_n1dn2', 'clone_src', 'clone_src_', 'clone_stats', ...
    'inds', 'trees', 'clone_intra_read', 'clone_mutations');
%%
save([tree_dir 'fastas.mat'], 'clone_fasta');



%% post process for silent/non_silent  --  OBSELETE

% find all the clones that are in frame
% junk = [ mod(clone_stats(:,3),3) [clone_mat.RF]' ];
% ix = find( junk(:,1) == 0 & junk(:,2) == 0)';
[~, I, J, K] = silent_map();

for i=1:N
    if mod(i,1000)==0, fprintf('%d ', i); end 
    if clone_stats(:,11) ~= 0, continue; end
    v = clone_ids(i,2);

    
    sequences = map(cell2mat({clone_fasta{i}.Sequence}'));
    consensus = mode(sequences(2:end,:),1);
    root = sequences(1,:);
    sequences = [root; consensus];

    % trim until first codon base
    Vtrim = mod(clone_stats(i,3),3);
    % convert to codons and trim from the right
    sequences   = sequences(:, Vtrim+1:end);
    sequences  = seqs2codons(sequences(:,1:end-mod(size(sequences,2),3)));

    
    % fill in junctions of germline with root of tree 
    sequences(1,sequences(1,:) == 65) = sequences(2,sequences(1,:) == 65);
    sequences(2,sequences(2,:) == 65) = sequences(1,sequences(2,:) == 65);
    sequences(sequences == 65) = 1;
    assert(max(sequences(:)) < 65);
    x = (sequences(1,:)-1)*64 + sequences(2,:);
    clone_mutations(i).single_base = sum(J(x(:))); % single-base mutation
    clone_mutations(i).silent =      sum(I(x(:))); % single-base silent mutation
    clone_mutations(i).third_base =  sum(K(x(:))); % single-base third base mutation    
    
    % show silent/nonsilent mutations normalized by lengh of whole sequence
    len = size(sequences,2);
    bins = ceil(30* (1:len)/ len );
    F = false(len,30);
    F(sub2ind(size(F), 1:len, bins)) = true;

    y = sum(I(x),1);    % silent
    clone_mutations(i).germ_silentspots = int8(y*F);        
    y = sum(J(x)-I(x),1);  % non-silent
    clone_mutations(i).germ_nonsilentspots= int8(y*F);   
    
end
fprintf('\n');
