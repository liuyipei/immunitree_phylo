% refind FR1 reads.
load([tree_dir 'fastas.mat']);
load('~/scratch/data/lymph/t5.fa.mat');
data = correct_PCR_artifacts(data);

%% Find the Ds of each read
rep = load_repertoire('ihmmune_collapsed');
[junk2, junk] = max(clone_src, [], 2);
[S R] = ind2sub([4 6], junk);
[~,d] = ismember(clone_n1dn2(:,3), {rep.D.Header});

map = []; map('ACGTNR') = [1:5 5];
dict = 'ACGTN';

ix = find(R<=2);  % only FR1 reads
%%
err = [];
b = struct('lastV', [], 'firstJ', [], 'in1', [], 'in2', [], 'eaten', [], 'germline', []);
b = b(ones(N,1),1);

for i=ix' % for each variant
    if mod(i,100) == 0, fprintf('%d\n',i); end
    [~,~,id] =  play_with_clone(clone_fasta{i});
    v_ = clone_ids(i,2); j_ = clone_ids(i,3);
    V_seq = rep.V(v_).Sequence;
    J_seq = rep.J(j_).Sequence;

    % align reads
    edge_V = 40; edge_J = 20;
    [L_aligned R R_aligned trimmed_VJ] = VJ_align(data.reads(id), V_seq, J_seq, [edge_V edge_J]);

    if length(R) > 1
        if length(unique(cellfun(@length, R))) > 1, err = [err i]; continue; end
    end
    reads = int16(cell2mat(R));        
    clone = [L_aligned reads R_aligned];

%    consensus = mode(double(clone),1);
    consensus = mode(double(reads),1);

    % construct the germline                    
    noise = 0.001;
    d_ = max(1,d(i)); % no d=0
    [~, ~, al] = profileHMM_align({consensus}, {V_seq(end-edge_V+1:end)}, {rep.D(d_).Sequence}, {J_seq(1:edge_J)}, noise);
    germline = map([V_seq(trimmed_VJ(1)+1:end-edge_V) dict(al.germline) J_seq(edge_J+1:end-trimmed_VJ(2))]);

        
    % indexes of junctions
    a = [];
    a.lastV  = length(V_seq)-trimmed_VJ(1)-al.eaten(1);
    a.firstJ = length(J_seq)-trimmed_VJ(2)-al.eaten(4);
    a.firstJ = length(germline)+1 - a.firstJ; % correct to index from left
    a.in1 = a.lastV  + ([1:length(al.N1)]);
    a.in2 = a.firstJ - ([length(al.N2):-1:1]); 
    a.eaten = [trimmed_VJ(1) al.eaten trimmed_VJ(2)];
    a.germline = germline;
    b(i) = a;

    new_fasta = clone_fasta{i};
    new_seqs = dict([germline; clone]);
    for j=1:size(new_seqs,1)
        clone_fasta{i}(j).Sequence = new_seqs(j,:);
    end
    clone_stats(i,1) = size(clone,1);
    clone_stats(i,2) = length(germline);
    clone_stats(i,3:8) = a.eaten;
    clone_stats(i,9:10) = [a.lastV a.firstJ];
    clone_stats(i,11) = mod(clone_stats(i,2) + clone_stats(i,3) + clone_stats(i,8) - 1, 3);
end


