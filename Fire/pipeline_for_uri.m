function [a h chain res] = pipeline_for_uri(fasta_file, germline_file, chain_type, save_outputs)
% Chain Types:
% 0: Either Kappa or Lambda
% 1: Heavy, VDJ
% 2: Lambda, VJ
% 3: Kappa, VJ
% 4: TCR alpha, VJ
% 5: TCR beta, VDJ
% 6: TCR gamma, VJ
% 7: TCR delta, VDJ

close all hidden

X = fastaread(fasta_file);
fprintf('Starting with %d reads.\n', length(X));

has_D_segment = ismember(chain_type, [1,5,7])
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

if isempty(germline_file)    

    % find the V-J assignment for every read:
    for k=1:length(X),
        if mod(k,100) == 0
            k
        end
        [v(k),j(k),err(k)] = analyze_VDJ(X(k).Sequence, rep);  % for 3 arguments, it does not care for the D
    end
    fprintf('Throwing away %d reads with alignment issues.\n', sum(err));
    X = X(~err);  % work only with reads without issues
    v = v(~err);
    j = j(~err);

    % force all reads to align to the consensus V,J
    % report "problem" if the consensus is less than 90% of reads.
    v_ = mode(v)  
    j_ = mode(j)
    if sum(v == v_) < 0.9*length(X), fprintf('The V segments of the reads have concensus under 0.9.\n'); end
    if sum(j == j_) < 0.9*length(X), fprintf('The J segments of the reads have concensus under 0.9.\n'); end
    rep.V = rep.V(v_);
    rep.J = rep.J(j_);

else % assume that it is given...
    % assume the V,J combination is given to us in a file
    % create a mock repertoire with only those V and J
    G = fastaread(germline_file);
    rep.V = G(1); rep.J = G(2);
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
force_same_cluster = true;

% get reads from that length
iy = find(Z==l);
reads = int16(cell2mat(R(iy)));
if force_same_cluster
    c_ = ones(size(reads,1),1);    
else    
    noise = 0.001;

    [uniqueR, unique_map, dup_map] = unique(reads,'rows');

    % find probability for each read to be generated from germline VDJ. 
    [profileHMM_score, d] = profileHMM_align(R(iy(unique_map)), ...
        {V_seq(end-edge_VJ(1)+1:end)}, {rep.D.Sequence}, {J_seq(1:edge_VJ(2))}, noise);

    % cluster reads to clones
    if size(reads,1)>1
        nIter = 30;
        c_ = cluster_reads_to_clones_no_indels(uniqueR, profileHMM_score, dup_map, nIter, noise);
    else 
        c_ = 1;
    end
end

% get most prominant clusters
M = hist(c_, 1:max(c_))
fprintf('%d cluster.  Max cluster has %d reads.\n', length(M), max(M));

% Get the clone with the most reads
[~, max_clone] = max(M)
m = max_clone;
iz = find(c_ == m);            

% construct the full reads in the clone
clone = [L_aligned(iy(iz),:) reads(iz,:) R_aligned(iy(iz),:)];

% construct the germline for the clone
consensus = mode(double(reads(iz,:)),1);
germline = extrapolate_germline(consensus, edge_VJ, trimmed_VJ, rep);

%-%% %%%%%%%%%%%%%%%%%        Build Tree     %%%%%%%%%%%%%%%%%%%%

    % set no. of iterations
    if size(clone,1) > 10, 
        nIter = 2000; % large trees get more iterations                    
    else
        nIter = 200;        
    end 
    nIter=120
    
    % set priors
    priors = struct('is_codon', false, 'is_decay', false, 'pi', germline.seq);                   
    priors.NT = generate_sticky_prior(1e-3, 4, 630);
    priors.R = generate_sticky_prior(2e-3, 4, 6300000); 
    priors.nClasses = 1;

    % set initial values (and also fixed values) for mutation/noise model
    init = [];
    init.R = [generate_sticky_prior(2e-3, 4, 1); ones(1,4)];                    
    init.mut_model.NT = generate_sticky_prior(1e-3, 4, 1);
                    
    % disable these moves in the MCMC chain
    control = struct('mut_params', false, 'sites', false, 'read_params', false);
    
    % run chain
    [~, chain,stats,stats2,ll] = infer_tree(clone, nIter, priors, [], [], init, control);
    report_traces(chain(end),stats2,stats)
    
    %  set a to be the highest scoring *canonic* tree (and
    %  overwrite previous 'best')
    [k, ~, best] = get_best_canonic_tree(priors, ll, chain);
    k
    
    % clean and make the tree prettier
    a = convert_phylo_tree_to_mutation_tree(best);
    a = add_germline_as_root(a, germline.seq);
    a = collapse_edges(a);
    a = order_nodes(a);
    
    % update the content of the junctions in the germline based on the tree
    dict = 'ACGTN';
    germline.n1 = dict(best.sequences(1,germline.in1));
    germline.n2 = dict(best.sequences(1,germline.in2));
    germline.header = [germline.header '_' germline.n1 '_' germline.n2];
    
%%  figures and outputs  
    str = '';
    
    %pre 0110
    h = visualize_tree(a.tree, a.sequences, a.t, clone, 1, germline);
    if iy(iz(1)) == 1 %&& strcmp(X(1).Header, 'germline')
        str = sprintf(' target is in node %d', a.t(1));
    end
    h.Nodes(1).Label = [h.Nodes(1).Label str];

    show_read_alignment_to_germline(clone, germline, [], a.t);    
    set(gcf, 'Position', [6 281 1126 246]);
    res.best = best; res.stats = stats; res.stats2 = stats2; res.ll = ll;
    res.a = a; res.rep = rep; res.germline = germline; res.priors = priors;
    
    res.best = best; res.stats = stats; res.stats2 = stats2; res.ll = ll;
    res.a = a; res.rep = rep; res.germline = germline; res.priors = priors;    
    if save_outputs
        [~, username] = system('whoami');
        liuyipei = isequal(username(1:8), 'liuyipei');
        
        jpg_file_name = [fasta_file '.1.jpg'];
        if liuyipei, jpg_file_name = regexprep(jpg_file_name, '.*/', './output/'); end
        saveas(gcf, jpg_file_name);   

        jpg_file_name = [fasta_file '.2.jpg'];   
        if liuyipei, jpg_file_name = regexprep(jpg_file_name, '.*/', './output/'); end
        saveas(h.hgAxes, jpg_file_name); 
        
        % save tree to file
        clone_file = [fasta_file '.mat'];
        if liuyipei, clone_file = regexprep(clone_file, '.*/', './output/'); end

        save(clone_file, 'clone', 'best', 'stats', 'stats2', 'll', 'a', 'rep', 'germline');
        fprintf('Saved results to files.\n');
    end

%%%%%%%%%%%%%%%%%%%%%        Done With Tree     %%%%%%%%%%%%%%%%%%%%




end
