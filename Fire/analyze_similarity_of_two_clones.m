function analyze_similarity_of_two_clones(clone1, clone2)    
    Z1 = fastaread(clone1);
    Z2 = fastaread(clone2)
    consensus1 = char(mode(double(cell2mat({Z1(2:end).Sequence}'))));
    consensus2 = char(mode(double(cell2mat({Z2(2:end).Sequence}'))));
    ix = find(Z1(1).Header == '_');
    junc_length = length(Z1(1).Header)-ix(end-1)-1;    
    germ1 = Z1(1).Sequence;
    germ1(germ1 == 'N') = Z1(1).Header([ix(end-1)+1:ix(end)-1 ix(end)+1:end]);
            
    fprintf('%s (reads: %d to-germ: %d)\n', clone1, length(Z1)-1, sum(germ1~=consensus1));
    fprintf('%s (reads: %d to-germ: %d)\n', clone2, length(Z2)-1, sum(germ1~=consensus2));

    fprintf('junc-length: %d  dist-between-consesnus of two clones: %d \n\n', junc_length, sum(consensus1 ~= consensus2));
    consensus = [consensus1; consensus2];
    mut =  consensus ~= [germ1; germ1];
    consensus(~mut) = '-'

    reads = [cell2mat({Z1(2:end).Sequence}'); cell2mat({Z2(2:end).Sequence}')];
    B = pdist(double(reads), 'hamming')*length(germ1);
    figure; imagesc(squareform(B))
    plot_class_lines([1 length(Z1) size(reads,1)], false, 'm3', true);
    plot_class_lines([1 length(Z1) size(reads,1)], true, 'm3', true);
    
    clone1(clone1 == '_') = '-'; clone2(clone2 == '_') = '-';
    title( sprintf([clone1(find(clone1=='/', 1, 'last')+1:end) '\n' clone2(find(clone2=='/', 1, 'last')+1:end)]));

end

%%  
function test()
close all
tree_dir = '/afs/cs/u/joni/scratch/data/lymph/t5.fa_results/';
fid = fopen('/afs/cs/u/joni/scratch/data/lymph/junctions_shared13_34.txt');
junk = textscan(fid, '%s %*s %*s %*s %*s %*s');
fclose(fid);
files = junk{1};
N = length(files);
j = 29;
while j<N        
    analyze_similarity_of_two_clones([tree_dir files{j}], [tree_dir files{j+1}]);       
    j = j+2;    
end

%%  Show read similarity



%%  convert my clone file to original reads
fasta_file_name = '/afs/cs/u/joni/scratch/data/lymph/t5.fa';
tree_dir = '/afs/cs/u/joni/scratch/data/lymph/t5.fa_results/';
clones = {'patient_3_V_202_J_10_len_79_clone_1.fa', 'patient_4_V_202_J_10_len_79_clone_1.fa', ...
    'patient_3_V_77_J_6_len_86_clone_2.fa', 'patient_4_V_77_J_6_len_86_clone_6.fa'};
clones = cellfun(@(x) [tree_dir x], clones, 'uniformoutput', false);
find_original_reads_of_clone(fasta_file_name, clones, data, big_fasta);


end

%% Show pairs of clones from the same individual that have the same variant region
function test2()
%%
close all
files = {'patient_10_V_229_J_6_len_86_clone_1.fa', 'patient_10_V_229_J_6_len_86_clone_2.fa', ...
         'patient_8_V_175_J_9_len_86_clone_1.fa', 'patient_8_V_175_J_9_len_86_clone_2.fa', ...
         'patient_7_V_217_J_6_len_89_clone_1.fa', 'patient_7_V_217_J_6_len_89_clone_4.fa', ...
         'patient_8_V_77_J_6_len_62_clone_1.fa', 'patient_8_V_77_J_6_len_62_clone_2.fa', ...
         'patient_9_V_60_J_13_len_82_clone_1.fa', 'patient_9_V_60_J_13_len_82_clone_2.fa', ...
         'patient_10_V_32_J_12_len_83_clone_2.fa', 'patient_10_V_33_J_12_len_83_clone_4.fa'};
N = length(files);
j = 1;
while j<N        
    analyze_similarity_of_two_clones([tree_dir files{j}], [tree_dir files{j+1}]);       
    j = j+2;    
end

end
