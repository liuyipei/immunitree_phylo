% parallelize VDJ
clear;
cd ~/JVL/src/phylo/Fire
addpath('.');
addpath('../../phylo/');
addpath('../../DPtrees/');
addpath('../../VDJ/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));
fasta_file_name = '/afs/cs/u/joni/scratch/data/lymph/t5.fa';
load([fasta_file_name '.mat'], 'data');

%% compute the V-D-J of every read

rep_code = 'ihmmune_collapsed';
map('ACGTN') = 1:5;
dict = 'ACGTN';
N = length(data.reads);
progress = ceil(N/100);
for j=1:100
    vdj = zeros(N,3, 'int16');
    a = struct('N1', [], 'N2', [], 'eaten', []);
    a = a(ones(N,1));
    temp_file = sprintf([fasta_file_name '.vdj%d'], j);
    if exist([temp_file '.start'], 'file'), continue; end
    system(['touch ' temp_file '.start']);
    for i=((j-1)*progress+1):min(j*progress,N)
        x = map(data.reads{i});
        try
            [vdj(i,1),vdj(i,3),al, err, vdj(i,2)] = analyze_VDJ(x, rep_code);
            a(i).N1 = dict(al.N1); 
            a(i).N2 = dict(al.N2);
            a(i).eaten = al.eaten;
        end
    end
    %data.vdj = vdj;
    save([temp_file '.mat'], 'vdj', 'a');
end
exit
%% regroup
N = 384462;
progress = ceil(N/100);
vdj = zeros(N,3, 'int16');
a = struct('N1', [], 'N2', [], 'eaten', []);
a = a(ones(N,1));
fprintf('[          ]\n[');
for j=1:100
    if mod(j,10)==0, fprintf('*'); end
    temp_file = sprintf([fasta_file_name '.vdj%d'], j);
    loaded = load([temp_file '.mat'], 'vdj', 'a');
    vdj = vdj + loaded.vdj;
    a(((j-1)*progress+1):min(j*progress,N)) = loaded.a(((j-1)*progress+1):min(j*progress,N));
end
fprintf('\n\n');



%%  for each read, get the alignment based on the  ihmmune VDJs
fasta_file_name = '/afs/cs/u/joni/scratch/data/lymph/t5.fa';
load([fasta_file_name '.mat']);      
%%  for each read, get alignments
rep = load_repertoire('ihmmune_collapsed');
N = length(data.reads);
trimmed_VJ = zeros(N,4);
trimmed_VJ(:,3) =  1;
trimmed_VJ(:,4) = -1;

st = zeros(4,N);

for i=1:N %find(data.ihmmune_collapsed(:,1)==10)' %
    if mod(i,10000) == 0, fprintf('%d ', i); end
    j = data.ihmmune_collapsed(i,3);
    v = data.ihmmune_collapsed(i,1);
    if v == 0 || j==0, continue; end
    edge_V = 40; edge_J = 20;
    try
        [L_aligned R R_aligned trimmed_VJ(i,1:2) st(:,i) trimmed_VJ(i,4)] = ...
            VJ_align(data.reads(i), rep.V(v).Sequence, rep.J(j).Sequence, [edge_V edge_J]);
        trimmed_VJ(i,3) = length([L_aligned R{1} R_aligned]);    
    catch
        trimmed_VJ(i,4) = 10;
    end
end
fprintf('Done.\n');

%%
save('trimmed_VJ_for_reads', 'trimmed_VJ', 'st');


%%

for i=1:1000, 
    if mod(i,100)==0, fprintf('%d\n', i); end; 
    j = data.ihmmune_collapsed(i,3);
    v = data.ihmmune_collapsed(i,1);    
    [v d j al err] = analyze_VDJ(data.reads{i}, rep, data.ihmmune_collapsed(i,1), data.ihmmune_collapsed(i,3)); 
    if err>0, 
        fprintf('i=%d err=%d\n', i, err); 
    end; 
end



%%  how many reads contribute to more than one V-J combination
junk = data.ihmmune_collapsed(:,[1 4:7]);
N = length(data.reads);
junk2 = false(N,58);
junk2(sub2ind(size(junk2), (1:N)', junk(:,1)+1)) = true;
junk2(sub2ind(size(junk2), (1:N)', junk(:,2)+1)) = true;
junk2(sub2ind(size(junk2), (1:N)', junk(:,3)+1)) = true;
junk2(sub2ind(size(junk2), (1:N)', junk(:,4)+1)) = true;
junk2(sub2ind(size(junk2), (1:N)', junk(:,5)+1)) = true;
junk2(:,1) = false;  % reset the '0s'
junk4 = sum(junk2,2);
sum(junk4 >1)
sum(junk4 >0)
junk2(junk4 == 0, :) = true;

junk = double(data.vdj_collapsed(:,1));
junk3 = false(N,58);
junk3(sub2ind(size(junk3), (1:N)', junk(:,1)+1)) = true;
ix = junk3 > junk2;
sum(ix(:)) % the number of reads where there is a disagreement 4128.

