ones(10)*ones(10); % might help with Matlab bug 961964 on some machines (http://stackoverflow.com/questions/19268293/matlab-error-cannot-open-with-static-tls)
%%  Add to path
clear
set(0,'DefaultTextFontname', 'Lucida Console')

% Set variable 'phylo_path' the absolute location of the immunitree_phylo 
% directory taken from git hub
% We then add all related directories as well.
phylo_path = '/scratch/liuyipei/gitrepo/immunitree_phylo/';   % this line needs to be manually set
fprintf('phylo path: %s\n', phylo_path)
cd([phylo_path 'Fire'])
addpath(phylo_path, '-end');
addpath([phylo_path 'Fire'], '-end'); % fixes issues with affinegapmex/nwalign
addpath([phylo_path 'VDJ'], '-end');
addpath([phylo_path 'util'], '-end');
addpath([phylo_path 'overwrite_matlab'], '-begin');

%% matlabpool open local 4
%%  Just the heavy
dbstop if error;
% provide absolute filenames to input in the variable 'files'.
input_dir = [phylo_path, 'sample_immunitree_data/']; % this line needs to be manually set
files = dir([input_dir '*.fasta'])
opts = struct(...
    'upgma_done_on_unique_reads', true, ...
    'display_figures',            false, ...
    'save_displayed_figures',     false, ...
    'trim_to_full_multialigned',  true, ...
    'trim_to_any_multialigned',   false);
for i = 1:length(files)
    filename = files(i).name
    chain_type = 1 % assume heavy chain input
    if ~isempty(strfind(filename, 'light')) % check if the phrase "light" is part of the file name
        chain_type = 0 % assume light chain input
    end    
    save_results = true
    pipeline_stage = 0 % can be changed to numbers greater than 0 to skip initial steps
    germline_file = [] % generally considered unknown
    nIter = 300 % number of MCMC iterations -- 100+ encouraged
    [a h chain] = pipeline_for_hiv(... # despite the name, this is not necessarily just for HIV
        [input_dir filename], ...
        germline_file, ...
        chain_type, save_results, ...
        pipeline_stage, opts, ...
        nIter);
end

% generate figures and text files indicating the structures of the trees
sample_pipeline_for_figures

%% all output are located in:
output_location = [phylo_path, 'Fire/output']