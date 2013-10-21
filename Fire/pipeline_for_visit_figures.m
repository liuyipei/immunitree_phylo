%%  Add to path
clear

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

files = [dir('./output/*.fa.mat') dir('./output/*.fasta.mat')] % generate figures from all mat files
[~, username] = system('whoami');
liuyipei = isequal(username(1:8), 'liuyipei');
close all hidden
close all
cutoff_datenum = datenum(2012, 7, 24, 1, 0, 0);
%dict = 'ACGTN'; % not the usual case
dict = 'ACGT-';


% make a legend figure first
legend_figure_bg = biograph([0 1; 0 0]);
h=view(legend_figure_bg);
h.Nodes(1).UserData.Distribution = 0;
h.Nodes(2).UserData.Distribution = [1 1 1 1 1 4; ...
                                    0 1 1 1 1 4; ...
                                    0 1 4 0 3 6; ...
                                    0 0 2 0 0 4; ...
                                    0 0 0 0 0 2];
set(h.Nodes(2), 'Label', ...
    sprintf([...
    'Numbers on edges denote count of mutations. \n',...
    'Nodes are labeled with IDs and their sizes. \n',...
    'Meaning of outer wheel slices, shown ccw, \nfrom dark blue to orange: \n', ...
    '  Y2-V1 Y2-V2 Y2-V3 Y3-V1 Y3-V2 Y3-V3\n' ...
    'Meaning of inner wheel Y3-V3 (orange) subslices ccw,\n from medium blue to orange: \n', ...
    '  IgM IgD IgG IgA IgE\n' ...
    ]));
h.CustomNodeDrawFcn = @(node) dist_node_draw(node);
jpg_file_name = './output/slices.legend.png';
f = get(h.hgAxes, 'Parent');
set(f, 'HandleVisibility', 'on');
set(f, 'Position', [1 1 800 600]);
f = get(h.hgAxes, 'Parent');
print(f, '-dpng', jpg_file_name);

%find(arrayfun(@(x)isequal(files(x).name, 'IGHV3-74.IGHJ3.10.BFI-0000391.fa.clusterline.0001.fa.mat'), 1:363))
for j=1:length(files)
    if files(j).datenum < cutoff_datenum
        fprintf('skipping %s (due to datenum cutoff)\n', files(j).name);
        continue
    end

    %% setup
    close all
    close all hidden
    load(['./output/' files(j).name])
    
    a_tree_parents = a.tree(:, 1);
    a_tree_depths = zeros(1, length(a_tree_parents));
    for node = 1:length(a_tree_parents)
        parent = a_tree_parents(node);
        if parent == 0            
            a_tree_depths(node) = 0;
        else
            a_tree_depths(node) = 1 + a_tree_depths(parent);
        end
    end    
    [widest_depth, widest_depth_width] = mode(a_tree_depths);
    [~, widest_sibling_class] = mode(a_tree_parents);
    % descriptive statistics: max(depth), max(out degree), 
    % mean(node depth), mean node depth weighted by size
    description_textfile = [fasta_file '.description.txt']
    description_textfile = regexprep(description_textfile, '.*/', './output/');        
    [description_file_handle description_file_msg]= fopen(description_textfile, 'w+');
    fprintf(description_file_handle, '%s, %d\n', 'Max depth', max(a_tree_depths));
    fprintf(description_file_handle, '%s, %d\n', 'Size of largest sibling class', widest_sibling_class);
    fprintf(description_file_handle, '%s, %.3f\n', 'Mean depth', mean(a_tree_depths));
    fprintf(description_file_handle, '%s, %.3f\n', 'Read-count weighted mean depth', a_tree_depths*a.tree(:, 2)/a.tree(1,3));
    fprintf(description_file_handle, '%s, %d\n', 'Number of reads', a.tree(1,3));

    fclose(description_file_handle);
    
    
    header_textfile = [fasta_file '.header.txt']
    header_textfile = regexprep(header_textfile, '.*/', './output/');
    headers = {X(iz).Header}';
    sequences = {X(iz).Sequence}';
    [header_file_handle header_file_msg]= fopen(header_textfile, 'w+');
    fprintf(header_file_handle, 'header, node, sequence\n');
    for i = 1:length(headers)
        fprintf(header_file_handle, '%s, %d, %s\n', headers{i}, int16(a.t(i)), sequences{i});
    end
    fclose(header_file_handle);
    
    % compute the index of the final V base in the references (reference V,
    % and reference VJ) such that both references are populated at that
    % index. this is used to annotate the root germline node
    map('ACGTNRacgtn-') = [1:5 5 1:5 5];
    germline_vj_and_v_int = int16(cell2mat(cellfun(@(x) map(x(:)'), ...
        {upgma_alignment((end-2):end).Sequence}', ... % remove the concatenated germline_vj and the v-sequence
        'UniformOutput', 0)));
    last_full_refVandVJalign_pos = find(~max(germline_vj_and_v_int==5,[],1), 1, 'last');
    V_num_bases_ignored = sum(... % find the number of non-gaps from prefix that was dropped
        ismember(upgma_germline(1:(first_full_multialign_pos-1)), 'ATCGatcg'));
    V_num_bases_left = sum(... % find the number of non-gaps from prefix that was dropped
        ismember(upgma_germline(first_full_multialign_pos:end), 'ATCGatcg'));
    V_trimmed_germline.mutcount_end_indx = last_full_refVandVJalign_pos;
    
    % determine labels, based on visit number and year
    per_read_labels = 8 * ones(length(X), 1); % default value is 8 (dark dark red for 'unknown')
    
    assert(length(a.t) == length(X), 'a.t does not have the same length as X!')
    timestamp_strings = {'Y2-V1','Y2-V2','Y2-V3','Y3-V1','Y3-V2','Y3-V3'};
    for vi=1:size(X,1)
        for vj = 1:length(timestamp_strings)                    
            if strfind(X(vi).Header, timestamp_strings{vj}),
                per_read_labels(vi) = vj;
                continue
            end
        end
    end
    
    
    class_strings = {'IgM', 'IgD', 'IgG', 'IgA', 'IgE'};
    
    % per isotype, per time point, per NODE
    per_read_sublabels = zeros(5, length(timestamp_strings), size(a.tree,1));     
    
    % per_read_labels: (3 (or 5?) isotypes, 6 timestamps) % just stick with 5?
    %                       Y2-V1 ... ... Y3-V3
    %   IgM             [                        ]
    %   IgD             [                        ]
    %   IgG             [                        ]
    %   IgA             [                        ]
    %   IgE             [                        ]
    class_map = [1 2 3 4 5]; % staying with 5
    for vi=1:size(X,1)
        curr_time = per_read_labels(vi);
        for vj = 1:length(class_strings)                    
            if strfind(X(vi).Header, class_strings{vj}),
                per_read_sublabels(class_map(vj), curr_time, a.t(vi)) = ...
                    per_read_sublabels(class_map(vj), curr_time, a.t(vi)) + 1;
                continue
            end
        end
    end
    h = visualize_tree(a.tree, a.sequences, a.t, ...
        V_trimmed_clone, 1, V_trimmed_germline, ...
        per_read_labels, per_read_sublabels);
    
    set(h.Nodes(1), 'Label', ...
        sprintf(['[1]\nGermline: %s\n',...
        '%d leading germline V bases were not used. ', ...
        '%d VJ bases were used in phylogeny.\n', ...
        'Germline V mutation count cutoff read index:%d\n%s'], ...
        concat_germline_vdj.Header, V_num_bases_ignored, V_num_bases_left, ...
        V_trimmed_germline.mutcount_end_indx, regexprep(fasta_file, '^.*/', '')));


    jpg_file_name = [fasta_file '.2.png'];
    if liuyipei, jpg_file_name = regexprep(jpg_file_name, '.*/', './output/'); end
    f = get(h.hgAxes, 'Parent');
    set(f, 'HandleVisibility', 'on');
    set(f, 'Visible', 'off');
    set(f, 'Position', [1 1 max(800, widest_depth_width * 70) max(600, (1+max(a_tree_depths))*80 )]);
    print(f, '-dpng', jpg_file_name);
    
    % save images
    show_read_alignment_to_germline(V_trimmed_clone, V_trimmed_germline, [], a.t);
    set(gcf, 'Position', [50 50 1000 1000]); 
    
    jpg_file_name = [fasta_file '.1.png'];
    if liuyipei, jpg_file_name = regexprep(jpg_file_name, '.*/', './output/'); end
    print(gcf, '-dpng', jpg_file_name);    

        
    node_textfile = [fasta_file '.node.txt']
    node_textfile = regexprep(node_textfile, '.*/', './output/');        
    [node_file_handle node_file_msg]= fopen(node_textfile, 'w+');
    fprintf(node_file_handle, 'node, parent, node_size, subtree_size, description, sequence\n');
    
    convert_padded_numeric_to_nuc = @(x) dict(x(x<5));
    for i = 1:length(h.Nodes)
        fprintf(node_file_handle, '%d, %d, %d, %d, %s, %s\n', i, a.tree(i,1), ...
            a.tree(i,2), a.tree(i,3), h.Nodes(i).Description, ...
            convert_padded_numeric_to_nuc(a.sequences(i,:)));
    end
    fclose(node_file_handle);
    
end
