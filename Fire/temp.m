%%
clear
cd /vision/u/liuyipei/joni/phylo/Fire
phylo_path = '/vision/u/liuyipei/joni/phylo/';
addpath(phylo_path, '-end');
addpath([phylo_path 'Fire'], '-end'); % fixes issues with affinegapmex/nwalign
addpath([phylo_path 'VDJ'], '-end');
addpath([phylo_path 'util'], '-end');
%% matlabpool open local 4

clear
dbstop if error;
%hiv_dir = '/vision/u/liuyipei/vdj_data/Uri0806fasta/';
hiv_dir = '/vision/u/liuyipei/vdj_data/boyd_t13t9_cluster94_fa/';
files = [dir([hiv_dir '*.fasta']) dir([hiv_dir '*.fa'])]

opts = struct('upgma_done_on_unique_reads', true, ...
              'display_figures',            true, ...
              'save_displayed_figures',     true, ...
              'trim_to_full_multialigned',  true, ...
              'trim_to_any_multialigned',   false); 
% false/false on both trim_to options would mean no trimming at all

%for i = 4:length(files)
for i = 1:3
    filename = files(i).name
    currsize = files(i).bytes
    if currsize > 6e6 % 6MB
        continue
    end
    chain_type = 1; % 1 = heavy
    if ~isempty(strfind(filename, 'light'))
        chain_type = 0; % light
    end
    chain_type

    save_outputs = true;
    pipeline_stage = 0;
    pipeline_for_hiv([hiv_dir filename], [], chain_type, save_outputs, ...
        pipeline_stage, opts);
end
