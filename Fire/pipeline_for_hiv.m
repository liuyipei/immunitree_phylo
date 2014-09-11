function [a h chain res] = pipeline_for_hiv(fasta_file, ...
  germline_file, chain_type, save_outputs, pipeline_stage, opts, nIter)
% Chain Types:
% 0: Either Kappa or Lambda
% 1: Heavy, VDJ
% 2: Lambda, VJ
% 3: Kappa, VJ
% 4: TCR alpha, VJ
% 5: TCR beta, VDJ
% 6: TCR gamma, VJ
% 7: TCR delta, VDJ

if nargin < 5
    pipeline_stage = 0;
    % we may wish to instead load the parsed fasta sequences
end
if nargin < 6
    opts = struct(...
        'upgma_done_on_unique_reads', true, ...
        'display_figures',            false, ...
        'save_displayed_figures',     false, ...
        'trim_to_full_multialigned',  true, ...
        'trim_to_any_multialigned',   false);
end

dict = 'ACGTN';
if pipeline_stage == 0
    X = fastaread(fasta_file);
    [uniq_reads, raw_from_uniq, uniq_from_raw] = uniqueRowsCA({X.Sequence}');
    
    fprintf('Starting with %d reads. (%d uniq reads)\n', ...
        length(X), length(raw_from_uniq));

    has_D_segment = ismember(chain_type, [1,5,7]);
    if isempty(germline_file)    
        
        if chain_type == 0
            rep = load_repertoire('rep/IGlight', false, 'VJ');
            rep.D.Sequence = '';
            rep.D.Header = 'NoD';
        elseif chain_type == 1
            rep = load_repertoire('rep/IGH');
        elseif chain_type == 2 
            rep = load_repertoire('rep/IGL', false, 'VJ');
            rep.D.Sequence = '';
            rep.D.Header = 'NoD';
        elseif chain_type == 3
            rep = load_repertoire('rep/IGK', false, 'VJ');
            rep.D.Sequence = '';
            rep.D.Header = 'NoD';
        end


        % find the V-J assignment for every read:
        % for k=1:length(X),
        v   = -ones(1, length(X)); v_per_ur   = -ones(1, length(raw_from_uniq));
        j   = -ones(1, length(X)); j_per_ur   = -ones(1, length(raw_from_uniq));
        err = -ones(1, length(X)); err_per_ur = -ones(1, length(raw_from_uniq));
        parfor ur=1:length(raw_from_uniq) % foreach uniq read
            [v_per_ur(ur), j_per_ur(ur), err_per_ur(ur)] = ...
                analyze_VDJ(uniq_reads{ur}, rep);  % analyze_VDJ on 3 outputs: ignores D
        end
        for k = 1:length(X) 
            % propagate the computed values to reads identical to the ones analyzed by analyze_VDJ
            % ur is the uniq read index
            ur = uniq_from_raw(k);
            v(k) = v_per_ur(ur);
            j(k) = j_per_ur(ur);
            err(k) = err_per_ur(ur);
        end
        
        %%OLD LOGIC  
        %   for k = 1:length(X) 
        %   [v_old(k), j_old(k), err_old(k)] = ...
        %   analyze_VDJ(X(k).Sequence, rep);  % analyze_VDJ on 3 outputs: ignores D
        %end        
  
        fprintf('Throwing away %d reads with alignment issues and keeping the remaining %d reads\n', sum(err~=0), sum(err ==0));
        X = X(~err);  % work only with reads without issues
        v = v(~err);
        j = j(~err);
        [uniq_reads, raw_from_uniq, uniq_from_raw] = uniqueRowsCA({X.Sequence}');

        % force all reads to align to the consensus V,J
        % report "problem" if the consensus is less than 90% of reads.
        v_ = mode(v)  
        j_ = mode(j)
        if sum(v == v_) < 0.9*length(X), fprintf('The V segments of the reads have concensus under 0.9.\n'); end
        if sum(j == j_) < 0.9*length(X), fprintf('The J segments of the reads have concensus under 0.9.\n'); end
        rep.V = rep.V(v_);
        rep.J = rep.J(j_);

    else
        % assume that it is given...
        % assume the V,J combination is given to us in a file
        % create a mock repertoire with only those V and J
        junk = fastaread(germline_file);
        rep.V = junk(1); rep.J = junk(2);
        
        err = zeros(1,length(X));
        v = ones(1,length(X));
        j = ones(1,length(X));
        v_ =1;
        j_ =1;
    end

    fprintf('pipeline stage 0 (fasta read and VJ-ihmmualign) finished!\n');
    parsed_fasta_file = [regexprep(fasta_file, '.*/', './parsed_fasta/') '.parse_only.mat']
    save(parsed_fasta_file, 'X', 'err', 'germline_file', ...
      'has_D_segment', 'j', 'j_', 'rep', 'v', 'v_', 'dict', ...
      'uniq_reads', 'raw_from_uniq', 'uniq_from_raw', 'opts');
else
    parsed_fasta_file = [regexprep(fasta_file, '.*/', './parsed_fasta/') '.parse_only.mat']
    load(parsed_fasta_file)
    fprintf('loaded %d parsed reads from mat: pipeline stage %d\n', length(X), pipeline_stage)
end

fprintf('Aligning reads to V-J parts.\n');
V_seq = rep.V.Sequence;
J_seq = rep.J.Sequence;
edge_VJ = [40 25];
[L_aligned R R_aligned trimmed_VJ] = VJ_align({X.Sequence}, V_seq, J_seq, edge_VJ);

% get most prominant lengths
Z = cellfun(@length,R);
h = hist(Z, 0:1:max(Z));
h = h(2:end); % do not include length zero
[~, max_length] = max(h);
l = max_length;
sparse(h)
fprintf('%d bins of length.  Max_length = %d.\n', sum(h>0), l);

%%  Clustering reads (of same length) into clones
force_same_cluster = true

% get reads from that length
% iy = find(Z==l);
% reads = int16(cell2mat(R(iy)));
iy = 1:length(R);
reads = cellfun(@int16, R(iy), 'UniformOutput', false);

if force_same_cluster
    cluster_count_max = 1
else
    cluster_count_max = Inf
end

use_dirchlet = false
if use_dirchlet ==1 
    noise = 0.001;
    [unique_reads_sorted, unique_map, dup_map] = uniqueRowsCA_filler(reads,[],'first');

    % find probability for each read to be generated from germline VDJ. 
    [profileHMM_score, d] = profileHMM_align(R(iy(unique_map)), ...
        {V_seq(end-edge_VJ(1)+1:end)}, {rep.D.Sequence}, {J_seq(1:edge_VJ(2))}, noise);

    % cluster reads to clones
    cluster_iterations = 30;

    [c_ seqs aligned_reads] = cluster_reads_to_clones(reads, profileHMM_score, ...
      dup_map, cluster_iterations, noise, cluster_count_max);

    % get most prominant clusters
    M = hist(c_, 1:max(c_))
    fprintf('%d cluster.  Max cluster has %d reads.\n', length(M), max(M));

    % Get the clone with the most reads
    [~, max_clone] = max(M)
    m = max_clone;
    iz = find(c_ == m);            
    fprintf('Dominant clone has %d reads.', length(iz))

    % construct the full reads in the clone

    clone = [L_aligned(iy(iz),:) cell2mat(aligned_reads(iz)) R_aligned(iy(iz),:)];
    dirchlet_consensus = seqs{m};
    V_trimmed_germline = extrapolate_germline(dirchlet_consensus, edge_VJ, trimmed_VJ, rep);
else
    consensus_file = [regexprep(fasta_file, '.*/', './parsed_fasta/') '.consensus.mat']
    if pipeline_stage <= 1
        map('ACGTNRacgtn-') = [1:5 5 1:5 5];
        fprintf('Assuming one cluster and high quality reads\n');
        iz = 1:length(X);
        
        concat_germline_vdj = struct(...
            'Header', [regexprep(rep.V.Header, '\|.*$', '') ' | ' rep.J.Header], ...
            'Sequence', [rep.V.Sequence rep.J.Sequence]);
        
        % From Matlab Documentation:
        %   When Seqs is a cell or char array, SeqsMultiAligned is a 
        %   char array with the output alignment following the same order as the input.
        %%% old: 
        
        if opts.upgma_done_on_unique_reads    
            %transitionmatrix_5x5 = generate_sticky_with_indel_prior(2e-3, 5, 1, 0.05); % based on SHM, rather than read errors
            %scoringmatrix_4x4 = log(transitionmatrix_5x5(1:4,1:4));
            %gapopen = log(transitionmatrix_5x5(1,end));
            %gapextend = gapopen;            
            
            nuc44_full = nuc44;    
            rand4x4 = rand(4); 
            rand4x4_symm = rand4x4+rand4x4';
            noise_regularization = rand4x4_symm * 100 * eps;
            scoringmatrix = nuc44_full(1:4,1:4) + noise_regularization;
            
            % HIV settings
            %Y = multialign([uniq_reads; rep.V.Sequence; concat_germline_vdj.Sequence], ...
            %    'UseParallel', true, 'Weights', 'equal', 'Verbose', false, ...
            %    'TerminalGapAdjust', true, 'ExistingGapAdjust', false, ...
            %    'GapOpen', 6.9329, 'ExtendGap', 6.9329-100*eps, ...
            %    'ScoringMatrix', scoringmatrix);                 
            Y = multialign([uniq_reads; rep.V.Sequence; concat_germline_vdj.Sequence], ...
                'UseParallel', true, 'Weights', 'equal', 'Verbose', false, ...
                'TerminalGapAdjust', true, 'ExistingGapAdjust', true, 'GapOpen', 2*6.933);
            
            upgma_alignment = [X; rep.V; concat_germline_vdj]; % just transfer the headers, in order
            for i = 1:length(uniq_from_raw)
                upgma_alignment(i).Sequence = Y(uniq_from_raw(i),:);
            end
            upgma_alignment(end-1).Sequence = Y(end-1,:);
            upgma_alignment(end).Sequence = Y(end,:);
        else
            upgma_alignment = multialign([X; rep.V; concat_germline_vdj], 'UseParallel', true, 'TerminalGapAdjust', true);
        end
        
        upgma_consensus = seqconsensus(upgma_alignment); % no longer used, but useful to have
        
        % int16 is imporant so that you can add properly during inference
        clone = int16(cell2mat(cellfun(@(x) map(x(:)'), ... % just the raw reads from the fasta file
          {upgma_alignment(1:end-2).Sequence}', ... % concatenated germline_vj and the v-sequence removed
          'UniformOutput', 0)));  
        germline_vj_and_v_int = int16(cell2mat(cellfun(@(x) map(x(:)'), ...
          {upgma_alignment((end-2):end).Sequence}', ... % remove the concatenated germline_vj and the v-sequence
          'UniformOutput', 0)));
            
        % check for stop codons in the raw reads
        check_for_stop_codons = true;
        if check_for_stop_codons
            stop_codons = [map('TAG'); map('TAA'); map('TGA')];
            Y_trimming = X(raw_from_uniq); % create some extra fields, based on unique reads
            Y_trimming(1).nuc2trim = []; 
            Y_trimming(1).stopIndices = [];
            curr_unpadded = cell(length(Y_trimming),1);
            padded_v_seq_mapped = map(upgma_alignment(end-1).Sequence);
            for i = 1:length(Y_trimming) % during infer tree, seqs are in number-space   
                
                % identify the reading frame of the sequences
                [Y_trimming(i).nuc2trim, Y_trimming(i).Sequence] ...
                    = padded_nucseq_to_codon_position(clone(i,:), padded_v_seq_mapped);

                % skip the gaps by restricting < 5
                % get the locations of stop indices -- and trim that many
                % letters off the beginning of curr_unpadded{i}
                curr_unpadded{i} = clone(i, clone(i,:) < 5);
                Y_trimming(i).stopIndices = find_stop_codons_given_unpadded_inframe(...
                    stop_codons, curr_unpadded{i}((1 + Y_trimming(i).nuc2trim):end));
            end

            % display these stop codons as problems
            index_of_reads_with_stop_codons = find(~cellfun(@isempty, {Y_trimming.stopIndices}));
            idread_count = histc(uniq_from_raw, 1:max(uniq_from_raw));
            for i = index_of_reads_with_stop_codons(:)'
                X_indx = raw_from_uniq(i);
                codon_letters = dict(curr_unpadded{i}(Y_trimming(i).stopIndices +  Y_trimming(i).nuc2trim));
                fprintf('seq %d (with %d idreads): %s : %s:pos %d + %d.\n', ...
                          X_indx, idread_count(i), codon_letters, Y_trimming(i).Header, ...
                          Y_trimming(i).stopIndices(1), Y_trimming(i).nuc2trim)
            end
            fprintf('stop codons: %d/%d.\n',length(index_of_reads_with_stop_codons), length(Y_trimming))
        else
            Y_trimming = [];
        end
        fprintf('pipeline stage 1 (alignment and stop codon check) finished!\n');

        save(consensus_file, 'iz', 'upgma_alignment', 'upgma_consensus', ...
            'map', 'concat_germline_vdj', 'clone', 'Y_trimming',...
            'uniq_reads', 'raw_from_uniq', 'uniq_from_raw', 'opts');
    else
        fprintf('loading upgma consensus from mat file; pipeline stage was %d\n', pipeline_stage)
        load(consensus_file);
    end
    
    %trim the leading regions of the Vs so that FR1/FR2 all start populated
    if      opts.trim_to_full_multialigned
        first_full_multialign_pos = find(~max(clone==5,[],1),1,'first'); % position, before trimming
    elseif  opts.trim_to_any_multialigned
        first_full_multialign_pos = find(~min(clone==5,[],1),1,'first'); % position, before trimming
    else    % no trimming
        first_full_multialign_pos = 1; % no trimming -- i.e. start with the first base
    end
    
    V_trimmed_clone = clone; % trim away by setting the earlier base reads into gaps
    V_trimmed_clone(:,1:(first_full_multialign_pos-1)) = 5; % 5 as in gap
    upgma_germline = upgma_alignment(end).Sequence;    

    V_trimmed_germline = struct('seq', map(upgma_germline),  ...
        'Sequence', upgma_germline(first_full_multialign_pos:end), ...
        'in1', [], 'in2', [],...
        'header', concat_germline_vdj.Header, ...
        'frame_shift', 0);    
end
    
%-%% %%%%%%%%%%%%%%%%%        Build Tree     %%%%%%%%%%%%%%%%%%%%

    if ~exist('nIter', 'var')
        nIter = 2000
    end
    % set no. of iterations
    %if size(V_trimmed_clone, 1) > 10, 
    %    nIter = 2000; % large trees get more iterations                    
    %else
    %    nIter = 200;        
    %end 
    %nIter = 5000;
    
    
    % set priors
    priors = struct('is_codon', false, 'is_decay', false, 'pi', V_trimmed_germline.seq, 'NT', [], 'R', []);

    priors.NT = generate_sticky_with_indel_prior(1e-2, 5, 6.3e2, 0.05);
    priors.R = generate_sticky_with_indel_prior(2e-3, 5, 6.3e6, 0.65); % was 6.3e6
    priors.nClasses = 1;

    % set initial values (and also fixed values) for mutation/readerror model
    init = [];
    init.mut_model.NT = generate_sticky_with_indel_prior(1e-2, 5, 1, 0.05);
    init.R = generate_sticky_with_indel_prior(2e-3, 5, 1, 0.65);
                    
    % disable these moves in the MCMC chain
    control = struct('mut_params', false, 'sites', false, 'read_params', false);
    
    % ignore positions that are uniform -- for performance reasons only
    % note that this will offset all log-likelihoods by some constant
    [Vt_mode, Vt_mode_freq] = mode(double(V_trimmed_clone), 1);
    Vt_pos_diverse = Vt_mode_freq < size(clone, 1); % different base characters possible
    if sum(Vt_pos_diverse) < 2
        Vt_pos_diverse(1:2) = 1; % mark two spots as "diverse" just to make sure the resulting object is a matrix
    end
    fprintf('%d/%d base positions are diverse.\n', sum(Vt_pos_diverse), size(clone,2));
    Vt_clone_diverse = clone(:, Vt_pos_diverse);
    priors_Vt_diverse = priors; 
    priors_Vt_diverse.pi = V_trimmed_germline.seq(Vt_pos_diverse);
    % Note that these disabled MCMC moves may require a subtlely different
    % verbal interpretation, but the effects intuitively cancel each other out, mostly
    %   mut_params, sites, read_params
    
    % This MCMC chain is based on diverse positions only -- August 2012
    rseed = RandStream('mcg16807','Seed',0); % make sure results are reproducible
    RandStream.setDefaultStream(rseed);
    [~, chain, stats, stats2, ll] = infer_tree( ...
        Vt_clone_diverse, nIter, priors_Vt_diverse, [], [], init, control);     
    
    % make MCMC chain (old way, involving uniform base positions)
    %[~, chain, stats, stats2, ll] = infer_tree( ...
    %    V_trimmed_clone, nIter, priors, [], [], init, control); 
    
    %  find the highest scoring *canonic* tree
    %  re-add in the uniform  bases
    tic;
    [k, canonical_ll, best] = get_best_canonic_tree(priors_Vt_diverse, ll, chain, clone);
    best = read_uniform_bases(best, Vt_pos_diverse, Vt_mode);
    fprintf('Optimal canonical_ll: %f at iteration k = %d/%d.\n',max(sum(canonical_ll)), k, nIter);
    fprintf('Search for canonical_ll took %fs.\n', toc);
    
    % plot figures
    report_traces(chain(end),stats2,stats, sum(canonical_ll));      
    
    
    % clean, readd the uniform bases, add germline and make the tree prettier
    collapse_identical = true;
    a = convert_phylo_tree_to_mutation_tree(best, collapse_identical, clone);    
    
    a = add_germline_as_root(a, V_trimmed_germline.seq);
    a = order_nodes(a);
    a = collapse_edges(a);
    a = order_nodes(a);    
    % insert the removed uniform base positions 
    % back into the tree object (a)
    
    assert_identical_reads_map_to_same_cell(raw_from_uniq, uniq_from_raw, a.t);
    
    fprintf('final phylo tree (a) computed.\n');
    
    % update the content of the junctions in the germline based on the tree
    dict = 'ACGT-';
    V_trimmed_germline.n1 = dict(best.sequences(1,V_trimmed_germline.in1));
    V_trimmed_germline.n2 = dict(best.sequences(1,V_trimmed_germline.in2));
    V_trimmed_germline.header = [V_trimmed_germline.header '_' V_trimmed_germline.n1 '_' V_trimmed_germline.n2];
    fprintf('germline sequence and header documented.\n');
    
%%
    res.best = best; res.stats = stats; res.stats2 = stats2; res.ll = ll;
    res.a = a; res.rep = rep; res.germline = V_trimmed_germline; res.priors = priors;
    fprintf('res computed.\n');
    
%%  save mat file
    if save_outputs        
        % [~, username] = system('whoami');
        liuyipei = 1; %isequal(username(1:8),'liuyipei'); % default to this behavior

        % save tree to file
        clone_file = [fasta_file '.mat'];
        if liuyipei, clone_file = regexprep(clone_file, '.*/', './output/'); end

        save(clone_file, 'chain', 'clone', ...
            'V_trimmed_clone', 'concat_germline_vdj', 'V_trimmed_germline', ...
            'best', 'stats', 'stats2', 'll', 'a', 'first_full_multialign_pos', ...
            'rep', 'upgma_germline', 'upgma_alignment', 'save_outputs', 'res', ...
            'fasta_file', 'nIter', 'uniq_reads', 'raw_from_uniq', 'uniq_from_raw', ...
            'X', 'iz', 'opts',...
            'Vt_mode', 'Vt_mode_freq', 'Vt_pos_diverse', 'priors_Vt_diverse');
    
        fprintf('Saved simulation results to mat files.\n');

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
    end   
    %%%%%%%%%%%%%%%%%%%%%        Done With Tree     %%%%%%%%%%%%%%%%%%%%
    
    %%% figures-related section
    display_figures = opts.display_figures
    if display_figures
        %% compute/show figures    
        % [~, username] = system('whoami');
        liuyipei = 1; %isequal(username(1:8), 'liuyipei');
        
        show_read_alignment_to_germline(V_trimmed_clone, ...
            V_trimmed_germline, [], a.t);
        V_num_bases_ignored = sum(... % find the number of non-gaps from prefix that was dropped
            ismember(upgma_germline(1:(first_full_multialign_pos-1)), [1 2 3 4]));
        V_num_bases_left = sum(... % find the number of non-gaps from prefix that was dropped
            ismember(upgma_germline(first_full_multialign_pos:end), [1 2 3 4]));        
        
        % compute the index of the final V base in the references (reference V,
        % and reference VJ) such that both references are populated at that
        % index. this is used to annotate the root germline node
        germline_vj_and_v_int = int16(cell2mat(cellfun(@(x) map(x(:)'), ...
            {upgma_alignment((end-2):end).Sequence}', ... % remove the concatenated germline_vj and the v-sequence
            'UniformOutput', 0)));
        last_full_refVandVJalign_pos = find(~max(germline_vj_and_v_int==5,[],1), 1, 'last');
        V_trimmed_germline.mutcount_end_indx = last_full_refVandVJalign_pos;
        
        h = visualize_tree(a.tree, a.sequences, a.t, ...
            V_trimmed_clone, 1, V_trimmed_germline);     
        
        
        set(h.Nodes(1), 'Label', ...
            sprintf(['[1]\nGermline: %s\n%d leading germline V bases were not used.\n', ...
            '%d VJ bases were used in phylogeny.\n%s'], ...
            concat_germline_vdj.Header, V_num_bases_ignored, V_num_bases_left, ...
            regexprep(fasta_file, '^.*/', '')));
        
        set(gcf, 'Position', [10 10 1000 1000]); 

        % save images and node file
        save_displayed_figures = opts.save_displayed_figures
        if save_displayed_figures
            png_file_name = [fasta_file '.1.png'];
            if liuyipei, png_file_name = regexprep(png_file_name, '.*/', './output/'); end
            saveas(gcf, png_file_name);
            png_file_name = [fasta_file '.2.png'];   
            if liuyipei, png_file_name = regexprep(png_file_name, '.*/', './output/'); end
            saveas(h.hgAxes, png_file_name);
            
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
    end

end

