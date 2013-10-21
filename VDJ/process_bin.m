function process_bin(v,j,bin_dir) 
%%  Load repertoire of D's and the chosen V and J genes
    rep.V = fastaread('V.fa');
    rep.J = fastaread('J.fa');
    rep.D = fastaread('D.fa');

%%  Prepare data
real_data = 1;
if real_data == 1
%  Load real data
    v=6, j=4
    V_seq = rep.V(v).Sequence; 
    J_seq = rep.J(j).Sequence;
    Z = fastaread(sprintf('/afs/cs/u/joni/scratch/data/Uri/church/bin_%d_%d.fasta', v, j));
    
    X = {Z.Sequence};    
    N = length(X);
    for i=1:N
        ix = find(Z(i).Header==',', 3);
        assert(length(ix) == 3);
        run_id(i) = Z(i).Header(ix(end)+2)-'0';
    end
    clear Z
else
% Generate Clones
    % choose V and J randomly.
    v = ceil(length(rep.V)*rand);
    j = ceil(length(rep.J)*rand);
    V_seq = rep.V(v).Sequence;
    J_seq = rep.J(j).Sequence;
    
    map('ACGTNacgtn') = [1:5 1:5];
    nClones = 10; % how many clones to generate
    clones = cell(nClones,1);
    fprintf('Generating Clones...\n');
    for s=1:nClones
        % choose a D gene
        d = ceil(length(rep.D)*rand);
        D_seq = rep.D(d).Sequence;
        [clones{s} V N1 D N2 J] = profileHMM_gen(V_seq, D_seq, J_seq);
    end

% Generate reads from a DP tree rooted in each clone
    fprintf('Generating reads...\n');
    N = 2000;
    noise = 0.005; indel = 0.003;
%    [X c] = gen_reads(clones, N, noise, indel);
    [X c t seqs T] = gen_reads(clones, N, noise, indel, 0);
   
end    
%%  Trim the edges of the reads so they start and end in perfect alignment
    fprintf('Trimming reads to CDR region...\n');
    Y = cell(N,1);
    Y_nt = cell(N,1);
    cut_v = length(V_seq)-40;
    [trimmed_aligned R] = trim_and_fix_reads(X, V_seq, 0, cut_v);
    [L_aligned R] = trim_and_fix_reads(R, V_seq(cut_v:end), 0, 20);    
    [R_aligned R] = trim_and_fix_reads(R, J_seq, 1, 20);
    dict = 'ACGTN';
    for i=1:length(R)
        Y{i} = [L_aligned{i} R{i} R_aligned{i}];
        Y_nt{i} = dict(Y{i});
    end
    
    % index duplicate reads
    fprintf('Resolving duplicates in CDR region... ');
    [~, unique_map, dup_map] = unique(Y_nt);
    uniqueY = Y(unique_map);
    fprintf('They make up %.1f%% of the reads, so we have %d unique reads.\n', ...
        100*(1-length(uniqueY)/length(Y)), length(uniqueY));
%%  Map reads to profile HMM
    fprintf('Mapping the reads to profile HMMs...\n');
    tic;
    [profileHMM_score profileHMM_id] = profileHMM_align(uniqueY, {V_seq(cut_v:end)}, {rep.D.Sequence}, {J_seq});
    toc;    
    
    if real_data == 1 
        filename = sprintf('profileHMM_clone_%d_%d.mat', v, j);
        fprintf('profileHMM alignment scores saved to %s\n', filename);
        save(filename, ...
            'profileHMM_score', 'profileHMM_id', 'dup_map');
    end
            
%% Cluster + align reads
    fprintf('clustering...\n')
    [c_, clones_, aligned, edt] = cluster_reads_to_clones(Y, profileHMM_score, dup_map, 30);
    M = hist(c_, 1:max(c_));
    if real_data == 0,  confusion_matrix(c, c_), end
    
%%  Run DPTree
    Z = cell(N,1);
    T_ = cell(length(M),1);
    seqs_ = cell(length(M),1);
    t_ = cell(length(M),1);
    for i=1:length(R)
        Z{i} = [trimmed_aligned{i} aligned{i}];
    end      
    
    params = [];
    params.me = 0.5;
    params.alpha = 0.5;
    params.epsilon = 0.005;
    params.sigma = 0.025;
    for s = 1:length(M)
        fprintf('processing clone %d\n', s);
        [T_{s} seqs_{s} t_{s}]= DPtrees_inference_binary(cell2mat(Z(c_==s)), 50, params);
    end

%%
addpath('../DPtrees/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));

%%
save(sprintf('results_clone_%d_%d.mat', v, j), 'c_', 'aligned', 'edt', ...
    'T_', 'seqs_', 't_', 'run_id');
%%
visualize(c_, aligned, edt, T_, seqs_, t_, run_id);
%%
if real_data == 0, visualize(c_, aligned, edt, T, seqs); end

%% Time point stuff


end
