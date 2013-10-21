function find_original_reads_of_clone(fasta_file_name, clone_files, data, big_fasta)

%%  convert my clone file to original reads
%fasta_file_name = '/afs/cs/u/joni/scratch/data/lymph/t5.fa';
%tree_dir = '/afs/cs/u/joni/scratch/data/lymph/t5.fa_results/';
%clone_file = [tree_dir 'patient_3_V_202_J_10_len_79_clone_1.fa'];

if ~exist('big_fasta', 'var') || isempty(big_fasta)
    fprintf('loading...');
    big_fasta = fastaread(fasta_file_name);
    fprintf('done\n');
end
output_clone_files = cellfun(@(x) [x '.original'], clone_files, 'uniformoutput', false);
%%
for k=1:length(clone_files)
    Z = fastaread(clone_files{k});
    headers = cell2mat({Z(2:end).Header}');
    headers = headers(:,1:8);
    [TF J] = ismember(headers, data.read_id, 'rows');
    assert(all(TF));
    reads = big_fasta(J);
    fastawrite(output_clone_files{k}, reads);
end
%read_id = cellfun(@(x) x(18:25), {big_fasta.Header}, 'uniformoutput', true);
