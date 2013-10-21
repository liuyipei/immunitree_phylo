%%  Add to path
clear
set(0,'DefaultTextFontname', 'Lucida Console')
phylo_path = '/filesystem/u/liuyipei/joni/phylo/';
if strfind(pwd(), '/visionnfs/')
    phylo_path = regexprep(phylo_path, '/filesystem/', '/visionnfs/')
elseif strfind(pwd(), '/vision/')
    phylo_path = regexprep(phylo_path, '/filesystem/', '/vision/')
elseif strfind(pwd(), '/scail/')
    phylo_path = regexprep(phylo_path, '/filesystem/', '/scail/')
end
fprintf('phylo path: %s\n', phylo_path)
cd([phylo_path 'Fire'])
addpath(phylo_path, '-end');
addpath([phylo_path 'Fire'], '-end'); % fixes issues with affinegapmex/nwalign
addpath([phylo_path 'VDJ'], '-end');
addpath([phylo_path 'util'], '-end');
addpath([phylo_path 'overwrite_matlab'], '-begin');
%addpath(phylo_path, '-bgein');
%addpath([phylo_path 'Fire'], '-begin'); % fixes issues with affinegapmex/nwalign
%addpath([phylo_path 'VDJ'], '-begin');
%addpath([phylo_path 'util'], '-begin');

%%
%matlabpool open local 12
%%
matlabpool open local 4

%% Jessica Finn
clear; close force all;
dbstop if error;
hiv_dir = '/scail/u/liuyipei/vdj_data/sirota_apeltsin_1108/';

opts = struct('upgma_done_on_unique_reads', true, ...
              'display_figures',            true, ...
              'save_displayed_figures',     false, ...
              'trim_to_full_multialigned',  true, ...
              'trim_to_any_multialigned',   false);

files = [dir([hiv_dir '*.fa']), dir([hiv_dir '*.fasta'])]
%for i = 1:length(files)
%    filename = files(i).name
filename = '158_IGHV7-4-1_IGHJ4_CGPPPFSSGRYPPGYW.fa'
    chain_type = 1; % heavy 
    if ~isempty(strfind(filename, 'light'))
        chain_type = 0; % light
    end
    chain_type

    save_results = true
    pipeline_stage = 1
    [a h chain] = pipeline_for_hiv([hiv_dir filename], ...
        [], ...
        chain_type, save_results, pipeline_stage);
%end
%%