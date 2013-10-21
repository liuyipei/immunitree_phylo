%%  Add to path
clear
set(0,'DefaultTextFontname', 'Lucida Console')
phylo_path = '/vision/u/liuyipei/joni/phylo/';
if strfind(pwd(), 'visionnfs'), phylo_path = regexprep(phylo_path, '/vision/', '/visionnfs/'),end
if strfind(pwd(), 'scail'), phylo_path = regexprep(phylo_path, '/vision/', '/scail/'),end
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
matlabpool open local 12
%%
matlabpool open local 4
%% 075 is a short one
%% without a screen, 005 
%% several sequences w diff lengths, 021
%% singleton node, 020
%% temporal 106--> 98 (before only, after only)
clear
dbstop if error;
hiv_dir = '/afs/cs/u/joni/scratch/data/Uri/HIV/';
files = dir([hiv_dir '098.fasta']);
[a h chain] = pipeline_for_hiv([hiv_dir files(1).name], [], 1, false, 2);

%%  Just the heavy
clear
dbstop if error;
%hiv_dir = '/vision/u/liuyipei/vdj_data/highIdHIV/';
%hiv_dir = '/vision/u/liuyipei/vdj_data/HIVset/';
hiv_dir = '/vision/u/liuyipei/vdj_data/Uri0725_cluster_heavy_fa/';
files = dir([hiv_dir '*.fa'])
for i = 1:length(files)
    filename = files(i).name
    chain_type = 1; % heavy 
    if ~isempty(strfind(filename, 'light'))
        chain_type = 0; % light
    end
    chain_type
    
    save_results = true
    pipeline_stage = 1        
    [a h chain] = pipeline_for_hiv([hiv_dir filename], ...
        ['/vision/u/liuyipei/vdj_data/Uri0725fasta/' 'PGT145_heavy_root.fasta'], ...
        chain_type, save_results, pipeline_stage);
end
%%
%%  Just the light
clear
dbstop if error;
%hiv_dir = '/vision/u/liuyipei/vdj_data/highIdHIV/';
%hiv_dir = '/vision/u/liuyipei/vdj_data/HIVset/';
hiv_dir = '/vision/u/liuyipei/vdj_data/Uri0725_cluster_light_fa/';
files = dir([hiv_dir '*.fa'])
for i = 1:length(files)
    filename = files(i).name
    chain_type = 1; % light 
    if ~isempty(strfind(filename, 'light'))
        chain_type = 0; % light
    end
    chain_type
    
    save_results = true
    pipeline_stage = 0        
    [a h chain] = pipeline_for_hiv([hiv_dir filename], ...
        ['/vision/u/liuyipei/vdj_data/Uri0725fasta/' 'PGT145_light_root.fasta'], ...
        chain_type, save_results, pipeline_stage);
end
%%

%% all in one file: 23584_heavy.broad.high-ident.fasta'
clear
dbstop if error;
hiv_dir = '/vision/u/liuyipei/vdj_data/Uri0725fasta/';
files = dir([hiv_dir '*.fasta']);

% 23584_heavy.broad.high-ident.fasta
% 23584_light.broad.high-ident.fasta
% PGT145_heavy_root.fasta
% PGT145_predicted_precursors.fasta
% PGT145_light_root.fasta

chain_type = 1; % heavy 
% if length(findstr(filename, 'light')) > 0
%     chain_type = 0 % light
% end
save_outputs = true;
pipeline_stage = 0;
pipeline_for_hiv([hiv_dir '23584_heavy.broad.high-ident.fasta'], [hiv_dir 'PGT145_heavy_root.fasta'], chain_type, save_outputs, pipeline_stage);

%% 23584_heavy.broad.high-ident.fasta
clear
dbstop if error;
hiv_dir = '/vision/u/liuyipei/vdj_data/Uri0725fasta/';
files = dir([hiv_dir '*.fasta']);

chain_type = 0;
save_outputs = true;
pipeline_stage = 0;
pipeline_for_hiv([hiv_dir '23584_heavy.broad.high-ident.fasta'], [hiv_dir 'PGT145_heavy_root.fasta'], chain_type, save_outputs, pipeline_stage);

%%

%% 23584_light.broad.high-ident.fasta
clear
dbstop if error;
hiv_dir = '/vision/u/liuyipei/vdj_data/Uri0725fasta/';
files = dir([hiv_dir '*.fasta']);

save_outputs = true;
pipeline_stage = 0;
chain_type = 1;
pipeline_for_hiv([hiv_dir '23584_light.broad.high-ident.fasta'], [hiv_dir 'PGT145_light_root.fasta'], chain_type, save_outputs, pipeline_stage);


%%%%%%%%%%%%%%%%%%%% while in mid execution after checking for internal
%%%%%%%%%%%%%%%%%%%% stop codons
X_vdj = [X; concat_germline_vdj];
[unique_X_vdj, m_indx_X_vdj, n_indx_X_vdj] = unique({X_vdj.Sequence});
glocal_scores = zeros(length(unique_X_vdj), length(unique_X_vdj));

for i = 1:length(unique_X_vdj)
   parfor j = i:length(unique_X_vdj)
       glocal_scores(i,j) = nwalign(unique_X_vdj{i}, unique_X_vdj{j}, 'GLOCAL', true);       
   end
   for j = i:length(unique_X_vdj)
       glocal_scores(j,i) = glocal_scores(i,j);
   end
end
glocal_eig = eig(glocal_scores)
[V D] = eig(glocal_scores)


all_edge_graph = glocal_scores + glocal_scores';
all_edge_graph = all_edge_graph - diag(diag(all_edge_graph));

wrong_directed_graph = graphminspantree(sparse(-all_edge_graph), 'method', 'Kruskal');
dists = graphallshortestpaths(-wrong_directed_graph, 'Directed', false);

root_node = n_indx_X_vdj(end); % the last one was the artificial germline root --> where did it map to after unique?
[~, nodes_ordered_from_root] = sort(dists(root_node,:));
topo_val_of_node_indexed_by_node = 1:length(unique_X_vdj);
topo_val_of_node_indexed_by_node(nodes_ordered_from_root) = 1:length(unique_X_vdj);

correct_directed_graph = zeros(length(unique_X_vdj),length(unique_X_vdj));
[ii,jj] = find(wrong_directed_graph);
for k = 1:length(ii);
    i = ii(k); j = jj(k);
    val = wrong_directed_graph(i,j);
    if topo_val_of_node_indexed_by_node(i) > topo_val_of_node_indexed_by_node(j)
        % i is the parent of j
        correct_directed_graph(j,i) = val;
    else
        correct_directed_graph(i,j) = val;
    end
end

bg = biograph(correct_directed_graph)
view(bg)

figure(1); clf;
plot(log(abs(glocal_eig))); hold on;
plot(log(abs(glocal_eig)));
title('eigenvalues are mostly positive');
xlabel('ascending eigenvalues');
ylabel('log(abs(x))');

clustergram(glocal_scores)

inv_gs = inv(glocal_scores);
gs_sorted = sort(inv_gs(:));
gs_thresh = gs_sorted(floor(0.9*length(gs_sorted)));
inv_gs_thresholded = (inv_gs > gs_thresh).*inv_gs;
inv_gs_thresholded = inv_gs_thresholded - diag(diag(inv_gs_thresholded));
inv_gs_th_bg = biograph(inv_gs_thresholded);


%% Boyd 2012/04/08

clear
dbstop if error;
hiv_dir = '/vision/u/liuyipei/vdj_data/boyd/';
files = dir([hiv_dir '*.fasta']);
i = 3
filename = files(i).name
chain_type = 1; % heavy 
[a h chain] = pipeline_for_hiv([hiv_dir filename], [], 1, false, 2);
%% Boyd: Chen Wang 2012/06/01
clear
dbstop if error;
% hiv_dir = '/vision/u/liuyipei/vdj_data/chenwang_06_2012/';
hiv_dir = '/visionnfs/u/liuyipei/vdj_data/chenwang_06_2012/';
files = dir([hiv_dir '*.fasta']);
opts = struct('upgma_done_on_unique_reads', true, ...
              'display_figures',            true, ...
              'save_displayed_figures',     true, ...
              'trim_to_full_multialigned',  true, ... # some reads are FR1 and others are FR2 in the chen dataset
              'trim_to_any_multialigned',   false);
for i = 1:length(files)
    close force all
    filename = files(i).name
    chain_type = 1; % heavy 
    save_outputs = true;
    pipeline_stage = 0;
    [a h chain] = pipeline_for_hiv([hiv_dir filename], [], chain_type, save_outputs, pipeline_stage);
end


%%
chain_type = 1;
fasta_file = '/vision/u/liuyipei/vdj_data/Jake_VH_VJ_cluster_fa/IGHV1-IGHJ3.len17_0000379.fa.clusterline.0002.fa'
pipeline_for_hiv(fasta_file, [], chain_type, false, 0)

%%
y2v1_flag = cellfun(@(x)isempty(findstr('Y2-V1', x)), {X.Header})
y2v2_flag = cellfun(@(x)isempty(findstr('Y2-V2', x)), {X.Header})
y2v3_flag = cellfun(@(x)isempty(findstr('Y2-V3', x)), {X.Header})
y3v1_flag = cellfun(@(x)isempty(findstr('Y3-V1', x)), {X.Header})
y3v2_flag = cellfun(@(x)isempty(findstr('Y3-V2', x)), {X.Header})
y3v3_flag = cellfun(@(x)isempty(findstr('Y3-V3', x)), {X.Header})
%%

[ur reads_from_ur ur_from_reads] = unique(reads, 'rows');
LLs_ur = zeros(2, length(reads_in_v));  % unique reads
for read_indx = 1:length(ur_from_reads)
    ur_indx = ur_from_reads(read_indx)
    LLs_ur(:,ur_indx) = LLs_ur(:,ur_indx) + LLs(:,read_indx);
end

%%

[uniq_reads, raw_from_uniq, uniq_from_raw] = uniqueRowsCA({X.Sequence}');
Y = multialign([uniq_reads; rep.V.Sequence; concat_germline_vdj.Sequence], 'UseParallel', true, 'TerminalGapAdjust', true, 'Weights', 'equal');
upgma_alignment = [X; rep.V; concat_germline_vdj]; % just transfer the headers, in order
for i = 1:length(uniq_from_raw)
    upgma_alignment(i).Sequence = Y(uniq_from_raw(i),:);
end
upgma_alignment(end-1).Sequence = Y(end-1,:);
upgma_alignment(end).Sequence = Y(end,:);
%%%
        % find the V-J assignment for every read:
        % for k=1:length(X),
        for k=raw_from_uniq
            [v(k),j(k),err(k)] = analyze_VDJ(X(k).Sequence, rep);  % for 3 outputs, it does not care for the D
        end
        for k = 1:length(X)
            % k_ belongs in raw_from_uniq, hence it has alraedy been populated
            k_ = raw_from_uniq(uniq_from_raw(k)); 
            v(k) = v(k_);
            j(k) = j(k_);
            err(k) = err(k_);
        end
%%
[uniq_reads reads_from_ur ur_from_reads] = unique(reads, 'rows');

uniq_reads_t_ = zeros(1,length(raw_from_uniq));
for i = 1:length(uniq_from_raw)
    x = uniq_from_raw(i);
    if uniq_reads_t(x) == 0
        uniq_reads_t_(x) = t_(i);
    else
        assert(uniq_reads_t_(x) == t_(i), 'Error: identical reads were mapped to different cells!')
    end
end


%% J-065.IGHV1-69_IGHJ6_449.fasta
clear; close force all;
dbstop if error;
hiv_dir = '/vision/u/liuyipei/vdj_data/Uri0806fasta/';

opts = struct('upgma_done_on_unique_reads', true, ...
              'display_figures',            true, ...
              'save_displayed_figures',     true, ...
              'trim_to_full_multialigned',  false, ...
              'trim_to_any_multialigned',   true);

chain_type = 1;
save_outputs = true;
pipeline_stage = 2
pipeline_for_hiv([hiv_dir 'J-065.IGHV1-69_IGHJ6_449.fasta'], [], ...
    chain_type, save_outputs, pipeline_stage, opts);



%% Uri1015, 80% id so that all pgts are included
clear; close force all;
dbstop if error;
hiv_dir = '/visionnfs/u/liuyipei/vdj_data/Uri1015_dnaclust80fa/';

opts = struct('upgma_done_on_unique_reads', true, ...
              'display_figures',            false, ...
              'save_displayed_figures',     false, ...
              'trim_to_full_multialigned',  false, ...
              'trim_to_any_multialigned',   true);

files = [dir([hiv_dir '*.fa']), dir([hiv_dir '*.fasta'])]
for i = 1:length(files)
    filename = files(i).name
    chain_type = 1; % heavy 
    if ~isempty(strfind(filename, 'light'))
        chain_type = 0; % light
    end
    chain_type
    
    save_results = true
    pipeline_stage = 0       
    [a h chain] = pipeline_for_hiv([hiv_dir filename], ...
        [], ...
        chain_type, save_results, pipeline_stage);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% plos revision

clear; close force all;
dbstop if error;
hiv_dir = '/scail/u/liuyipei/vdj_data/Uri_PlosPath_revision/';

opts = struct('upgma_done_on_unique_reads', true, ...
              'display_figures',            false, ...
              'save_displayed_figures',     false, ...
              'trim_to_full_multialigned',  false, ...
              'trim_to_any_multialigned',   true);

files = [dir([hiv_dir '*.fa']), dir([hiv_dir '*.fasta'])]
for i = 1:length(files)
    filename = files(i).name
    chain_type = 1; % heavy 
    if ~isempty(strfind(filename, 'light'))
        chain_type = 0; % light
    end
    chain_type
    
    save_results = true
    pipeline_stage = 0       
    [a h chain] = pipeline_for_hiv_plospathrevision([hiv_dir filename], ...
        [], ...
        chain_type, save_results, pipeline_stage);
end