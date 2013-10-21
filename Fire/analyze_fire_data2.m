%%  Load parameters
clear;
cd ~/JVL/src/phylo/Fire
addpath('.');
addpath('..');
addpath('../VDJ/');
addpath('../util/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));

%%
fasta_file_name = '/afs/cs/u/joni/scratch/data/lymph/t5.fa';
load([fasta_file_name '.mat']);      

% random seed is initialized based on current time
defaultStream = RandStream.getDefaultStream;
load('random_seed.mat');
random_seed(1) = uint32(1e13*mod(now,1e-5));
defaultStream.State = random_seed;

%%  Approach #2
method = 'ihmmune_collapsed';
rep = load_repertoire(method, false); % true for ihmmune

if ~exist([fasta_file_name '_results'], 'dir')
    system(['mkdir ' fasta_file_name '_results']);
end

% DONE: preprocess data to account for "jackpots"
data = correct_PCR_artifacts(data);

map = []; map('ACGTNR') = [1:5 5];
dict = 'ACGTN';

%min_size_of_cluster = 2;
register_non_clones = true;
min_size_of_clone = 1;
build_tree_for_single_replica_cluster_as_well = false;

%% inject variants mutation model 
rep = load_repertoire('ihmmune_collapsed');
A = get_IGH_alignments(rep, 'V');
B = get_IGH_alignments(rep, 'J');
load('NT_variants');
uninformed_prior = true

%% for each individual
for k = randperm(10) % 1:10
%%
    patient_file = sprintf('%s_results/patient_%d', fasta_file_name, k);
    if exist([patient_file '.begin'], 'file'), continue; end
    system(['touch ' patient_file '.begin']);

    % for each sample
    loc = 1:4;
    cur_sample = mod((k-1)*4+loc -1, 40)+1

%%
    % collect most prominant v-j combinations (at least XX reads)        
    [v j] = collect_VJs(data, get_reads(data,cur_sample, true), method, min_size_of_clone);

    % for each v-j combination
    for vj = randperm(length(v)) % 1:length(v) % increase parallelisation.
        v_ = v(vj), j_ = j(vj), k
        rep_ = rep; rep_.V = rep.V(v_); rep_.J = rep.J(j_);

        % qsub data
        pause(2*rand);  % to decrease chance of concurrency
        VJ_file = sprintf('%s_results/patient_%d_V_%d_J_%d', fasta_file_name, k, v_, j_);
        if exist([VJ_file '.start'], 'file'), continue; end
        system(['touch ' VJ_file '.start']);
        if exist([VJ_file '.fa'], 'file'), system(['rm ' VJ_file '.fa']); end


        V_seq = rep.V(v_).Sequence;
        J_seq = rep.J(j_).Sequence;
        ix = get_reads(data, cur_sample, true, [], [], [v_ 0 j_], method);            
        assert(length(ix) >= min_size_of_clone);        
        labels =  [mod(data.subject_num(ix)-1,4)+1 2*(double(data.FR(ix)-1)) + data.rep(ix)-'a'+1];
        
        % align reads
        edge_VJ = [40 20];
        [L_aligned R R_aligned trimmed_VJ] = VJ_align(data.reads(ix), V_seq, J_seq, edge_VJ);

        % get most prominant lengths
        Z = cellfun(@length,R);
        h = hist(Z, 0:1:max(Z));
        h = h(2:end); % do not include length zero

        noise = 0.001; % important:  Chance of a mismatch *compared to gemline* due to hyper-mutation or sequencing noise.            
        lengths_with_some_mass = find(h>=min_size_of_clone); 
        if isempty(lengths_with_some_mass), continue; end

        % for each length (with enough reads)
        for l = lengths_with_some_mass

            VJL_file = sprintf('%s_len_%d', VJ_file, l);

            % get reads from that length
            iy = find(Z==l);

            reads = int16(cell2mat(R(iy)));
            [uniqueR, unique_map, dup_map] = unique(reads,'rows');

            % find probability for each read to be generated from
            % germline VDJ. TODO: if D is given, use that value. On
            % second thought - D is inaccurate.
            [profileHMM_score, d] = profileHMM_align(R(iy(unique_map)), ...
                {V_seq(end-edge_VJ(1)+1:end)}, {rep.D.Sequence}, {J_seq(1:edge_VJ(2))}, noise);

            % cluster reads to clones
            if size(reads,1)>1
                c_ = cluster_reads_to_clones_no_indels(uniqueR, profileHMM_score, dup_map, 30, noise);
            else 
                c_ = 1;
            end

            % get most prominant clusters
            M = hist(c_, 1:max(c_));

            %% for each clone (with enough reads) 
            total_clones = 0;
            for m=find(M>=min_size_of_clone)
                iz = find(c_ == m);            

                % construct the full reads in the clone
                clone = [L_aligned(iy(iz),:) reads(iz,:) R_aligned(iy(iz),:)];

                % collect data about the reads as strings
                read_idx = ix(iy(iz));  
                read_ids = mat2cell([read_idx double(labels(iy(iz),:))], ones(length(iz),1), 3);
                read_ids = cellfun(@(x) sprintf('%06d_source_%d_rep_%d', x), read_ids, 'UniformOutput', false);

                % construct the germline                    
                consensus = mode(double(reads(iz,:)),1);
                germline = extrapolate_germline(consensus, edge_VJ, trimmed_VJ, rep_)

                multi_replica_cluster = size(unique(labels(iy(iz),:), 'rows'), 1) > 1;

%%%%%%%%%%%%%%%%%%%%%        Build Tree     %%%%%%%%%%%%%%%%%%%%
                if multi_replica_cluster || build_tree_for_single_replica_cluster_as_well

                    total_clones = total_clones+1;
                    clone_file = sprintf('%s_clone_%d', VJL_file, total_clones);

                    % Compute tree over the reads in the clone
                    nIter = 200;
                    if size(clone,1) > 10, nIter = 2000; end % large trees get more iterations                    
                    priors = struct('is_codon', false, 'is_decay', false, 'pi', germline.seq);
                    
                    init = [];
                    priors.NT = generate_sticky_prior(1e-3, 4, 630);
                    priors.R = generate_sticky_prior(2e-3, 4, 6300000); % 1e-7?
                    init.R = [generate_sticky_prior(2e-3, 4, 1); ones(1,4)];
                    
                    if uninformed_prior
                        priors.nClasses = 1;
                        init.mut_model.NT = generate_sticky_prior(1e-3, 4, 1);
                    else                        
                        priors.nClasses = size(clone,2);
                        NTV = NT.V(:,:, A(v_,:));
                        NTV = NTV(:,:,eaten(1)+1 : end-eaten(2));

                        NTJ = NT.J(:,:, B(j_,:));
                        NTJ = NTJ(:,:,eaten(5)+1 : end-eaten(6));

                        init.mut_model.NT = generate_sticky_prior(1e-3, 4, 1);
                        init.mut_model.NT = init.mut_model.NT(:,:,ones(priors.nClasses,1));
                        init.mut_model.NT(:,:,1:lastV) = NTV;
                        init.mut_model.NT(:,:,firstJ:end) = NTJ;
                        init.mut_model.rate_class = 1:priors.nClasses;
                    end
                    
                    control = struct('mut_params', false, 'sites', false, 'read_params', false);

                    [~, chain,stats,stats2,ll] = infer_tree(clone, nIter, priors, [], [], init, control);

                    %  set a to be the highest scoring *canonic* tree (and
                    %  overwrite previous 'best')
                    [~, ~, best] = get_best_canonic_tree(priors, ll, chain);
                    a = convert_phylo_tree_to_mutation_tree(best);
                    a = add_germline_as_root(a, germline.seq);
                    a = collapse_edges(a);
                    a = order_nodes(a);

                    % save tree to file
                    save([clone_file '.mat'], 'best', 'stats', 'stats2', 'll', 'a', 'v_', 'j_', 'd', 'method', 'cur_sample', 'trimmed_VJ');
                    % uncomment to save all the chain
%                     chain = chain(1:100:end);
%                     save([clone_file '_all.mat'],  'chain');

                    % overwrite the default values read-cluster
                    cell_ids = a.t;
                    germline.n1 = dict(best.sequences(1,germline.in1));
                    germline.n2 = dict(best.sequences(1,germline.in2));

                    fasta_clone_file = [clone_file '.fa'];
                    if exist(fasta_clone_file, 'file'), system(['rm ' fasta_clone_file]); end
                else % not a clone
                    % default values
                    cell_ids = zeros(size(clone,1),1); % empty
                    fasta_clone_file = [VJ_file '.fa'];
                end
%%%%%%%%%%%%%%%%%%%%%        Done With Tree     %%%%%%%%%%%%%%%%%%%%

                % fill junctions using inferred values.  Also report
                % junctions starting indexes
                germline.header = [germline.header '_' germline.n1 '_' germline.n2];                

                % save read-cluster to file. (append if not a clone)
                for i=1:length(read_idx)
                    read_ids{i} = [data.read_id(read_idx(i),:) '_' read_ids{i} '_cell_' sprintf('%03d', cell_ids(i))];                        
                end  

                read_ids = [ {germline.header}    ; read_ids];
                fastawrite(fasta_clone_file, read_ids, dict([germline.seq; clone]));                      

            end                
        end
        system(['touch ' VJ_file '.done']);
    end
    system(['touch ' patient_file '.end']);
end
%%
exit
