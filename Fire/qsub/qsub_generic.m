%%  Set path and working dir 
%      python script needs to prepend value for fasta_dir and filename
phylo_path = '/filesystem/u/liuyipei/joni/phylo/';
if strfind(pwd(), '/vision/'), phylo_path = regexprep(phylo_path, '/filesystem/', '/vision/'),end
if strfind(pwd(), '/visionnfs/'), phylo_path = regexprep(phylo_path, '/filesystem/', '/visionnfs/'),end
if strfind(pwd(), '/scail/'), phylo_path = regexprep(phylo_path, '/filesystem/', '/scail/'),end
fprintf('phylo path: %s\n', phylo_path)
cd([phylo_path 'Fire'])

addpath(phylo_path, '-end');
addpath([phylo_path 'Fire'], '-end'); % fixes issues with affinegapmex/nwalign
addpath([phylo_path 'VDJ'], '-end');
addpath([phylo_path 'util'], '-end');
%%
%files = [dir([fasta_dir '*.fasta']) dir([fasta_dir '*.fa'])]; % fasta files and fa files
%file = files(i)  % i definition is prepended by python script
%filename = file.name
chain_type = 1; % heavy 
if length(findstr(filename, 'light')) > 0 % use the file name to determine whether light or heavy
    chain_type = 0 % light
end
pipeline_stage = 0
opts = struct('upgma_done_on_unique_reads', true, ...
              'display_figures',            true, ...
              'save_displayed_figures',     true, ...
              'trim_to_full_multialigned',  true, ...
              'trim_to_any_multialigned',   false);
[a h chain] = pipeline_for_hiv([fasta_dir filename], [], chain_type, true, pipeline_stage, opts);
exit
    
