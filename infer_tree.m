%function [best chain stats stats2 ll] = infer_tree(reads, nIter, priors, filename, labels, init, control)
% run the tree algorithm.
% Inputs: 
%  reads: a set of equal-length reads, given as an integer matrix.
%  nIter: no. of iterations (default 7500). 
%  priors: a structure holding the parameters for the probabilitic model, 
%          and also the germline sequence that will be used for the root of 
%          the tree.  If missing or empty, use default value.  Also, any
%          missing field in the structure will use the default value.
%  filename: if not empty, intermediate results are written to this file.
%  labels: ground-truth cluster assignments.  If given, algorithm reports
%          the RAND index as a score.
%  init: a structure holding initial values for the state variables.
%        Missing values are initialized using `initialize_chain'.
%  control: a structure that states which MCMC moves are "on" when we run 
%           the chain.
%
% Outputs:
%   best: the state with the highest likelihood
%   cur : the last state
%   stats and stats2 are some statistics taken during the run of the chain.
%   They can be given to the function 'report_traces'.
%    stats2(:,j) = [log_likelihood, total_mutations, total_mismatches,...
%                   total_cells, birth_rate, death_rate, decay_factor, ...
%                   time_since_start, rand_score];
%    stats(i,:,j) = [total_cells live_cells total_tree_length];
%   ll: shows the log-likelihood.
% ll = [ sum log R(x_j , y_{s_j}) ;      % read-to-seq 
%        -M * log(N_alive)        ;      % read-assignment-to-cell 
%        BD-tree-likelihood       ;
%        sum log M(y_i , y_{s_i}) ;      % seq-to-parent-seq
%        prior-over-mutation-model;
%        site-assignment-to-class ;
%       ];
% 
function [best chain stats stats2 ll] = infer_tree(reads, nIter, priors, filename, labels, init, control)

if ~exist('control_', 'var'), control_ = []; end
if ~exist('init', 'var'), init = []; end
if ~exist('labels', 'var'), labels = []; end
if ~exist('filename', 'var'), filename = []; end
if ~exist('nIter', 'var'), nIter = 7500; end    

%%% random_seed: uncomment to use the same random seed.
% defaultStream = RandStream.getDefaultStream;
% load('random_seed.mat');
% defaultStream.State = random_seed;

% set priors
if ~exist('priors', 'var'), priors = []; end
priors = init_priors(priors);

% initialize chain
leaf_level = 0;
F_ = 50; 
[tree_ sequences_ t_ F_ rates_ mut_model_ R_ ] = ...
    initialize_chain(reads, leaf_level, priors, F_, init);

% State variables:
%  F_: the snapshot time. Cells alive at this time are considered "alive".
%  tree_: each row represents a single cell in the tree, and states its
%         parent node, birth time, and death time. If the cell is alive,
%         its death_time = F_.
%  sequences_: matrix of all cell sequences, ordered by their cell index.
%  t_: assignments of reads to cells in the tree.  
%  rates_: a structure with birth and death rates for tree.
%  mut_model: a structure with all the parameters governing the mutatio
%             model.
%  R_: a 5x4 transition matrix for the read noise model. the last row is
%      all ones, it's there to support reads with gaps.


is_codon = isfield(mut_model_, 'AA');  % are we working on codons or nts?

% set all MCMC moves to "on", and let the user override
control_ = struct('rates', true, 'manip', true, 'aux', true, ...
'birth_death', true, 'tune_death', true, 'seqinf', true, ...
'seqinf2', true, 'mut_params', true, 'sites', true, 'read_params', true, ...
'read_ass', true, 'read_gappings', false, ...
'collapse_nodes_with_idseq', false, 'greedy', false);
control = struct_override(control_, control);

% when greedy is -1, some of the internal functions change behavior
% in which they start taking only moves that increase the likelihood
global greedy; greedy = control.greedy;

% init statistics structures
ll = zeros(6,nIter);
stats = zeros(4,3, nIter);
stats2 = zeros(9, nIter);
best.ll = -1e1000; % practically -Inf

% stores all states in the chain
chain = struct('tree', tree_, 'sequences', sequences_, 't', t_, ...
    'mut_model', mut_model_, 'R', R_, 'rates', rates_, ...
    'F', F_, 'll', [], 'iteration', 0, 'control', control);

% construct the full transition matrix from the mutation model
Qs = get_Qs(mut_model_);


% remember which reads are identical, for assert/checks
[~, raw_from_uniq, uniq_from_raw] = unique(reads, 'rows');
read_partition = struct( ...
    'raw_from_uniq', raw_from_uniq, ...
    'uniq_from_raw', uniq_from_raw, ...
    'uniq_weights',  histc(uniq_from_raw, 1:length(raw_from_uniq)));

all_timer = tic;  % start the timer
fprintf('Running for %d iterations. \n', nIter);

for j=1:nIter
    annealing_prob = max(min(2*(j/nIter)^2-1+eps, 1), 0);

    assert_inference_invariants(raw_from_uniq, uniq_from_raw, t_, tree_, F_, sequences_);
    tic
    if control.birth_death % spawn a new node, delete existing node
        [tree_ sequences_ t_ codes] = MH_birth_death(tree_, sequences_, F_, Qs, rates_, t_, reads, R_);
        fprintf('--->  Birth-Death: iteration %d, %d proposals.  %.2f percent accepted.  %.2f seconds\n', ...
            j, sum(codes), 100*sum(codes(1:2))/sum(codes(1:3)), toc);     
        codes; tic; 
    end    
    stats(1,:,j) = [size(tree_,1) sum(tree_(:,3) <F_) sum(tree_(:,3) -tree_(:,2))];
    
    
    assert_inference_invariants(raw_from_uniq, uniq_from_raw, t_, tree_, F_, sequences_);
    if control.manip % local manipulation of tree
        tree_ = fix_births(tree_, F_);
        [tree_, ~, codes] = MH_manipulate(tree_, sequences_, Qs, rates_, t_); 
        fprintf('--->  Local Manip: iteration %d, %d proposals.  %.2f percent accepted.  %.2f seconds\n', ...
            j, sum(codes), 100*sum(codes(1:4))/sum(codes), toc); 
        codes; tic; 
    end
    stats(2,:,j) = [size(tree_,1) sum(tree_(:,3) < F_) sum(tree_(:,3) -tree_(:,2))];
    
    assert_inference_invariants(raw_from_uniq, uniq_from_raw, t_, tree_, F_, sequences_);
    if control.aux  % global manipulation of tree
        [tree_, ~, codes] = MH_auxiliary_graph(tree_, sequences_, Qs, t_); 
        fprintf('--->  Auxil Graph: iteration %d, %d proposals.  %.2f percent accepted.  %.2f seconds\n', ...
            j, sum(codes), 100*codes(1)/sum(codes), toc); 
        codes; tic; 
    end
    stats(3,:,j) = [size(tree_,1) sum(tree_(:,3) < F_) sum(tree_(:,3) -tree_(:,2))];
    
    assert_inference_invariants(raw_from_uniq, uniq_from_raw, t_, tree_, F_, sequences_);
    if control.tune_death % change death time of nodes
        [tree_, codes] = MH_tune_death(tree_, F_, rates_, t_); 
        fprintf('--->  Tune  Death: iteration %d, %d proposals.  %.2f percent accepted.  %.2f seconds\n', ...
            j, sum(codes), 100*sum(codes(1:3))/sum(codes), toc); 
        codes; tic; 
    end
    stats(4,:,j) = [size(tree_,1) sum(tree_(:,3) < F_) sum(tree_(:,3) -tree_(:,2))];

    assert_inference_invariants(raw_from_uniq, uniq_from_raw, t_, tree_, F_, sequences_);
    if control.seqinf % infer sequences for all tree cells based on reads
        % only some of the tree nodes have a descendant with a read
        % associated with it.  Sample the sequences for those nodes only.
        % The other nodes are not tied to any evidence. Set all their 
        % sequence letters to zero.
        annealed = rand(1) < annealing_prob;
        sequences_ = gibbs_sequences_of_nodes_with_evidence(tree_, ...
            sequences_, F_, Qs, mut_model_.pi, t_, R_, reads, annealed,...
            read_partition);
        fprintf('--->EvidSeqInfer1: iteration %d, %.2f seconds\n', j, toc); tic;
        
        % For each maximal subtree of nodes not tied to evidence, propose 
        % to replace the entire subtree with one sampled from the prior.  
        % If proposal accepted - great. if not, keep the original subtree
        % and sample only its sequences from the prior (always accepted).
        [tree_ sequences_ t_ codes] = MH_generate_subtree_from_node_with_no_evidence(...
            tree_, sequences_, t_, rates_, F_, Qs, mut_model_.pi);
        fprintf('--->NoevSeqInfer1: iteration %d, %.2f seconds\n', j, toc); tic;
    end
    
    tic;
    if control.mut_params  % samples the mutation model
        [mut_model_ codes total_mutations ll(4:5,j)] = ...
            MH_mutation_parameters(tree_, sequences_, t_, mut_model_, priors);
        fprintf('--->  mut params : iteration %d, %d proposals.  %.2f percent accepted.  %.2f seconds\n', ...
            j, sum(codes(:,2)), 100*sum(codes(:,1))/sum(codes(:,2)), toc); tic;

    else % just compute likelihood and stats
        [ll(4:5,j) total_mutations] = MH_mutation_parameters(tree_, ...
            sequences_, t_, mut_model_, priors);
    end

    tic;
    if control.sites % sample association of sites to mutation class
        [mut_model_.rate_class, ll(6,j)] = update_site_assignments(mut_model_, tree_, sequences_);  
        Qs = get_Qs(mut_model_); % recompute full transition matrix
        fprintf('---> update sites: iteration %d, %.2f seconds\n', j, toc); 
    end
    
    % prepare a nucleotide version of the sequences (if working in codon mode)
    if is_codon
        sequences_nt_ = codons2seqs(sequences_);    
    else
        sequences_nt_ = sequences_;    
    end
        
    if control.read_params % sample the read noise model
        [R_ total_mismatches] = estimate_read_noise_parameters(sequences_nt_, reads, t_, priors.R);
    else % just the stat
        [~, total_mismatches] = estimate_read_noise_parameters(sequences_nt_, reads, t_, priors.R);
    end
    
    assert_inference_invariants(raw_from_uniq, uniq_from_raw, t_, tree_, F_, sequences_);
    tic;
    if control.read_ass % sample read assignments to cells in the tree
        [ll(1:2,j) t_] = slice_read_assignments(tree_, sequences_nt_, ...
            t_, F_, R_, reads, read_partition, annealing_prob);            
        fprintf('--->  Read Assign: iteration %d, %.2f seconds\n', j, toc); tic;
        
    else % just likelihood
        ll(1:2,j) = slice_read_assignments(tree_, sequences_nt_, t_, F_, R_, reads);            
    end
    assert_inference_invariants(raw_from_uniq, uniq_from_raw, t_, tree_, F_, sequences_);
    
    if control.collapse_nodes_with_idseq
        % merge nodes having identical sequences; 
        % match the gappings in the moved reads to the destination node;
        
        [tree_, sequences_, t_, reads, collapse_status] = ...
            tree_collapse_identical_nodes(tree_, sequences_, t_, reads, F_);    
        fprintf('---> Node Merges: %d nodes merged; %d children moved\n', ...
            collapse_status.nodes_merged, collapse_status.children_moved);
    end
    
    %% do we adjust read gappings?
    if control.read_gappings
        reads = adjust_read_gappings(tree_, sequences_, t_, reads, raw_from_uniq, uniq_from_raw);
    end
    
    % Do a second iteration of the sequence inference
    % assert_identical_reads_map_to_same_cell(raw_from_uniq, uniq_from_raw, t_); % assert that identical reads belong together
    assert_inference_invariants(raw_from_uniq, uniq_from_raw, t_, tree_, F_, sequences_);
    
    if control.seqinf2 % infer sequences for all tree cells based on reads
        % only some of the tree nodes have a descendant with a read
        % associated with it.  Sample the sequences for those nodes only.
        % The other nodes are not tied to any evidence. Set all their 
        % sequence letters to zero.
        annealed = rand(1) < annealing_prob;
        sequences_ = gibbs_sequences_of_nodes_with_evidence(tree_, ...
            sequences_, F_, Qs, mut_model_.pi, t_, R_, reads, annealed,...
            read_partition);
        fprintf('--->EvidSeqInfer2: iteration %d, %.2f seconds\n', j, toc); tic;
        
        % For each maximal subtree of nodes not tied to evidence, propose 
        % to replace the entire subtree with one sampled from the prior.  
        % If proposal accepted - great. if not, keep the original subtree
        % and sample only its sequences from the prior (always accepted).
        [tree_ sequences_ t_ codes] = MH_generate_subtree_from_node_with_no_evidence(...
            tree_, sequences_, t_, rates_, F_, Qs, mut_model_.pi);
        fprintf('--->NoevSeqInfer2: iteration %d, %.2f seconds\n', j, toc); tic;
    end    
    assert_inference_invariants(raw_from_uniq, uniq_from_raw, t_, tree_, F_, sequences_);
    
    % sample birth and death rates
    [rates_ ll(3,j)] = MH_rates(tree_, F_, rates_, priors, control.rates);
	fprintf('Resampled Rates: BirthR %.4f, DeathR %.4f\n', rates_.birth, rates_.death);         
    
    % put the entire state in a struct
    a = struct('tree', tree_, 'sequences', sequences_, 't', t_, ...
        'mut_model', mut_model_, 'R', R_, 'rates', rates_, ...
        'F', F_, 'll', sum(ll(:,j)), 'iteration', j,...
        'control', control);

    % show rand index if ground truth is given
    if ~isempty(labels)
        a_ = convert_phylo_tree_to_mutation_tree(a);
        rand_score = confusion_matrix(labels, a_.t);
    else
        rand_score = 0;
    end
    
    % compute and print some stats
    fprintf('====================================================================================\n');
    fprintf('%s  Likelihood: %.2f  cells: %d (+%d) mutations: %d   mismatches:  %d  RAND=%.3f\n', ...
        filename, sum(ll(:,j)), size(tree_,1), size(tree_,1)-stats2(4,max(j-1,1)), ...
        total_mutations, total_mismatches, rand_score);
    fprintf('====================================================================================\n');
    stats2(:,j) = [sum(ll(:,j)) total_mutations, total_mismatches, size(tree_,1), ...
        rates_.birth, rates_.death, mut_model_.decay, toc(all_timer), rand_score];

    % if this state has the best likelihood so far, store it at 'best'.
    if sum(ll(:,j))>best.ll
        best = a;
    end    
    
    % add state to the chain
    chain(end+1) = a;

    % save state (if desired)
    if ~isempty(filename) && mod(j, 100)== 0
        %report_traces(a, stats2, stats);
        fprintf('saving state...');
        save([filename '.mat'], 'best', 'stats2', 'stats', 'chain', 'reads');
        fprintf('Done.\n');
    end    
end 

end

function assert_inference_invariants(raw_from_uniq, uniq_from_raw, t_, tree_, F_, sequences_)
    
    asserts_enabled = false; % disabled for performance reasons
    if asserts_enabled

        % assert that identical reads belong together
        assert_identical_reads_map_to_same_cell(raw_from_uniq, uniq_from_raw, t_); 

        % assert nodes with reads attached are alive
        assert(all(tree_(unique(t_),3)==F_));

        % assert that subclone sequences are populated with nucs or gaps
        assert(isempty(find(sequences_(:) == 0, 1)));
    end
end

function updateGreedyViaAnneling(j, nIter) % didn't work well (based on LL), -Yi, Aug 2012
global greedy;
greedy = (2*j/nIter - 1) > rand(1);
fprintf('greedy? : %d\n', greedy)
end

function test()
%%
    b = load('../synthetic_01_21_11')
    X = fastaread('/afs/cs/u/joni/JVL/src/phylo/tandy.fa');
    map('ACGTN') = 1:5;
    reads = int16(map(cell2mat({X.Sequence}')));
    
    nIter = 100;
    
    % set priors
    priors = struct('is_codon', false, 'is_decay', false);                   
    priors.NT = generate_sticky_prior(1e-3, 4, 630);
    priors.R = generate_sticky_prior(2e-3, 4, 630); 
    priors.nClasses = 1;
    
    % run chain
    [~, chain,stats,stats2,ll] = infer_tree(reads, nIter, priors, [], [], [], []);
    report_traces(chain(end),stats2,stats)
    
    %  set a to be the highest scoring *canonic* tree (and
    %  overwrite previous 'best')
    [k, ~, best] = get_best_canonic_tree(priors, ll, chain);
    k
        
    % clean and make the tree prettier
    a = convert_phylo_tree_to_mutation_tree(b.truth);
    a = collapse_edges(a);
    a = order_nodes(a);

    h = view_tree(a.tree, a.sequences, 1, [a.t b.labels], reads);
    for t=1:length(h.Nodes), set(h.Nodes(t),'Label',sprintf('%d',a.tree(t,2))); end
end
