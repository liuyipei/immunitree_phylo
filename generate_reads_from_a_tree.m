run_as_a_job = 1;

%% path
cd ~/JVL/src/phylo/
addpath('.');
addpath('../DPtrees/');
addpath('../VDJ/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));
global codon2aa codon2nt nt2codon; [codon2aa codon2nt nt2codon]= get_maps();

figure_positions = [0           0           0           0;
        1084         645         721         165;
        607         544        1021         953
        1815         645         681         165;
        1404         894        1088         603];
RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',sum(100*clock)));

dir = '~/scratch/data/Fire_FL'; 
%% starting point:  Load/generate reads
if exist('priors', 'var'), clear; end
close all;
nClasses = 3;
priors.R = generate_sticky_prior(0.0001, 4, 630000);
priors.NT = generate_sticky_prior(0.001, 4, 630);
priors.AA = generate_sticky_prior(1/21, 21, 210);
priors.decay = [0 1e-3];  % mean and std
priors.pi = 100*ones(1,64);

if ~exist('samples', 'var') && exist('truth', 'var') && isfield(truth,'samples')
    samples = truth.samples;
end

if exist('samples', 'var') && ~isempty(samples)
    [reads labels truth] = generate_data('real', samples);
elseif exist('samples','var') && isempty(samples)
    load('synthetic_02_23_11');
else
    %                      generate_data(  str,  nReads, L, maxCells, maxTime, nClasses, br,  dr)       
    [reads labels truth] = generate_data('synth', 1000, 33, 4200,     10,      nClasses, 0.8, 0.5, priors);
end
strSamples = sprintf('%d_', truth.samples);
    

%% visualize ground truth (erasing dead subtrees)
if isfield(truth, 'tree')
    a = convert_phylo_tree_to_mutation_tree(truth);
    visualize_tree(a.tree, codons2seqs(a.sequences), [a.t labels], reads);
end


%% init MCMC (all we need is 'reads').
job_id = ceil(rand*100000);
defaultStream = RandStream.getDefaultStream;
load('random_seed.mat');
defaultStream.State = random_seed;

if isfield(truth, 'pi') && size(truth.pi,1) == 1 
    priors.pi = truth.pi;
end

if ~exist('leaf_level', 'var'), leaf_level = 0; end
if ~exist('F_', 'var'), F_ = 50; end
[tree_ sequences_ t_ F_ rates_ mut_model_ R_ ] = ...
    initialize_chain(reads, nClasses, leaf_level, priors, F_);

% mut_model = best.mut_model;  tree_ = best.tree; sequences_ = best.sequences; t_ = best.t; R_ = best.R; F_ = best.F;  rates_ = best.rates;
%   rate_class_ = truth.rate_class;
%   tree_=truth.tree; sequences_ = truth.sequences; t_ = truth.t; F_ = truth.F; 
%   R_ = truth.R; 
%   rates_  = truth.rates;
%   AA_ = truth.AA; NT_ = truth.NT; pi_ = truth.pi;
%NT_ = generate_sticky_prior(0.00001, 4, 100000)/100000; AA_ = generate_sticky_prior(1/21, 21, 1);
%%  inference
global greedy; greedy = 0;
% control = struct('rates', false, 'manip', false, 'aux', false, ...
%     'birth_death', true, 'tune_death', false, 'seqinf', false, ...
%     'mut_params', false, 'sites', false, 'read_params', false, ...
%     'read_ass', false);
  control = struct('rates', true, 'manip', true, 'aux', true, ...
      'birth_death', true, 'tune_death', true, 'seqinf', true, ...
      'mut_params', true, 'sites', true, 'read_params', true, ...
      'read_ass', true);

if ~exist('nIter', 'var'),nIter = 7500; end
fprintf('Running for %d iterations. \n', nIter);
nInner = 1;    
ll = zeros(6,nIter);
%stats = zeros(4,nInner, nIter);
stats = zeros(4,3, nIter);
stats2 = zeros(9, nIter);

% store all intermediate states
cur = struct('tree', tree_, 'sequences', sequences_, 't', t_, ...
    'mut_model', mut_model_, 'R', R_, 'rates', rates_, ...
    'F', F_, 'll', [], 'iteration', 0);

best.ll = -1e10; % low number
Qs = get_Qs(mut_model_);
all_timer = tic;
j = 0;

for j=j+1:nIter
    tic
    if control.birth_death
        [tree_ sequences_ t_ codes] = MH_birth_death(tree_, sequences_, F_, Qs, rates_, t_, reads, R_);
        fprintf('--->  Birth-Death: iteration %d-%d, %d proposals.  %.2f percent accepted.  %.2f seconds\n', j,i, sum(codes), 100*sum(codes(1:2))/sum(codes(1:3)), toc);     
        codes, tic; %stats(1,i,j) = size(tree_,1); 
    end
    stats(1,:,j) = [size(tree_,1) sum(tree_(:,3) <F_) sum(tree_(:,3) -tree_(:,2))];
    assert(isempty(find(sequences_(:) == 0, 1)));

    if control.manip
        tree_ = fix_births(tree_, F_);
        [tree_, ~, codes] = MH_manipulate(tree_, sequences_, Qs, rates_, t_); 
        fprintf('--->  Local Manip: iteration %d-%d, %d proposals.  %.2f percent accepted.  %.2f seconds\n', j,i, sum(codes), 100*sum(codes(1:4))/sum(codes), toc); 
        codes, tic; 
    end
    stats(2,:,j) = [size(tree_,1) sum(tree_(:,3) < F_) sum(tree_(:,3) -tree_(:,2))];
    assert(isempty(find(sequences_(:) == 0, 1)));

    if control.aux
        [tree_, ~, codes] = MH_auxiliary_graph(tree_, sequences_, Qs, t_); 
        fprintf('--->  Auxil Graph: iteration %d-%d, %d proposals.  %.2f percent accepted.  %.2f seconds\n', j,i, sum(codes), 100*codes(1)/sum(codes), toc); 
        codes, tic; 
    end
    stats(3,:,j) = [size(tree_,1) sum(tree_(:,3) < F_) sum(tree_(:,3) -tree_(:,2))];
    assert(isempty(find(sequences_(:) == 0, 1)));

    if control.tune_death
        [tree_, codes] = MH_tune_death(tree_, F_, rates_, t_); 
        fprintf('--->  Tune  Death: iteration %d-%d, %d proposals.  %.2f percent accepted.  %.2f seconds\n', j,i, sum(codes), 100*sum(codes(1:3))/sum(codes), toc); 
        codes, tic; 
    end
    stats(4,:,j) = [size(tree_,1) sum(tree_(:,3) < F_) sum(tree_(:,3) -tree_(:,2))];

    if control.seqinf
        sequences_ = gibbs_sequences_of_nodes_with_evidence(tree_, sequences_, F_, Qs, mut_model_.pi, t_, R_, reads);
        toc; tic;
        % TODO: remove the pi parameter from the call to this funciton
        [tree_ sequences_ t_ codes] = MH_generate_subtree_from_node_with_no_evidence(tree_, sequences_, t_, rates_, F_, Qs, mut_model_.pi);
        fprintf('--->  Seqs Infere: iteration %d-%d, %.2f seconds\n', j,i, toc); 
        assert(isempty(find(sequences_(:) == 0, 1)));
    end
    
    tic;
    if control.mut_params
%        [AA_ NT_ decay_ pi_ codes total_mutations ll(4:5,j)] = ...
        [mut_model_ codes total_mutations ll(4:5,j)] = ...
            MH_mutation_parameters(tree_, sequences_, t_, mut_model_, priors);
        fprintf('--->  mut params : iteration %d, %d proposals.  %.2f percent accepted.  %.2f seconds\n', j, sum(codes(:,2)), 100*sum(codes(:,1))/sum(codes(:,2)), toc); 
    else
        [total_mutations ll(4:5,j)] = MH_mutation_parameters(tree_, sequences_, t_, mut_model_, priors);
    end

    if control.rates
        [rates_ ll(3,j)] = MH_rates(tree_, F_, rates_, control.rates);
    end
            
    tic;
    if control.sites
%        [rate_class_, ll(6,j)] = update_site_assignments(rate_class_, AA_, NT_, decay_, pi_, tree_, sequences_);  
        [mut_model_.rate_class, ll(6,j)] = update_site_assignments(mut_model_, tree_, sequences_);  
        Qs = get_Qs(mut_model_);
        fprintf('---> update sites: iteration %d, %.2f seconds\n', j, toc); 
    end
    
    sequences_nt_ = codons2seqs(sequences_);    
    if control.read_params
        [R_ total_mismatches] = estimate_read_noise_parameters(sequences_nt_, reads, t_, priors.R);
    else
        [~, total_mismatches] = estimate_read_noise_parameters(sequences_nt_, reads, t_, priors.R);
    end
    
    tic;
    if control.read_ass 
%           [t_ ll(1:2,j)] = gibbs_read_assignments(tree_, sequences_nt_, F_, R_, reads);                   
            [ll(1:2,j) t_] = slice_read_assignments(tree_, sequences_nt_, t_, F_, R_, reads);            
        fprintf('--->  Read Assign: iteration %d, %.2f seconds\n', j, toc); 
    else
        ll(1:2,j) = slice_read_assignments(tree_, sequences_nt_, t_, F_, R_, reads);            
    end
    
    a = struct('tree', tree_, 'sequences', sequences_, 't', t_, ...
        'mut_model', mut_model_, 'R', R_, 'rates', rates_, ...
        'F', F_, 'll', sum(ll(:,j)), 'iteration', j);

    % show rand index
    if isfield(truth, 'tree')
        a_ = convert_phylo_tree_to_mutation_tree(a);
        rand_score = confusion_matrix(labels, a_.t);
    else
        rand_score = 0;
    end
    
    fprintf('====================================================================================\n');
    fprintf('samples: %s  job_id: %d  Likelihood: %.2f  cells: %d (+%d) mutations: %d   mismatches:  %d  RAND=%.3f\n', ...
        strSamples, job_id, sum(ll(:,j)), size(tree_,1), size(tree_,1)-stats2(4,max(j-1,1)), total_mutations, total_mismatches, rand_score);
    fprintf('====================================================================================\n');
    stats2(:,j) = [sum(ll(:,j)) total_mutations, total_mismatches, size(tree_,1), rates_.birth, rates_.death, mut_model_.decay, toc(all_timer), rand_score];

    if sum(ll(:,j))>best.ll
        best = a;
    end    
    
    if mod(j, 100)== 0
        report_traces(a, stats2, stats);
        saveas(figure(5), sprintf('%s/runs/pictures/S%s%d.jpg', dir, strSamples, job_id));

%         figure(6); plot(stats2(8,1:j-1));
        cur(end+1) = a;
        fprintf('saving state...');
        save(sprintf('%s/runs/S%s%d.mat', dir, strSamples, job_id), 'best', 'stats2', 'stats', 'truth', 'cur', 'job_id', 'll');
        fprintf('Done.\n');
    end    
end 
truth
leaf_level
samples
if exist('run_as_a_job', 'var')
    if length(samples) == 1
	fid = fopen(sprintf('%s/jobs/S%d.done', dir, samples), 'w');
	fprintf(fid, sprintf('S%s%d\n', strSamples, job_id));
	fclose(fid);
    end
    exit
end

if 0
%% visualize best reconstruction (erasing dead subtrees)
a = best;
%%
a = best;
a = convert_phylo_tree_to_mutation_tree(a);
h = visualize_tree(a.tree, codons2seqs(a.sequences), [a.t labels], reads);
a = best;

%%
figure(4);
if isfield(truth, 't')                                  
    report_scores(truth, a, labels);
else
    report_scores(reads, a);
end
%%
report_traces(a, stats2, stats);


%%
saveas(h.hgAxes, sprintf('pictures/S%s%d_1.jpg', strSamples, job_id));
for k=[2 3 5 4]
    set(figure(k), 'position', figure_positions(k,:));
    saveas(figure(k), sprintf('pictures/S%s%d_%d.jpg', strSamples, job_id, k));
end

if isfield(truth, 't')
    a = convert_phylo_tree_to_mutation_tree(truth);
    h = visualize_tree(a.tree, codons2seqs(a.sequences), [a.t labels], reads);
    saveas(h.hgAxes, sprintf('pictures/S%s%d_0.jpg', strSamples, job_id));
    a = best;
end
    
%%  Show particles
for i=1:5:length(cur);
    a=cur(i);
    a_ = convert_phylo_tree_to_mutation_tree(a);
    h = visualize_tree(a_.tree, codons2seqs(a_.sequences), [a_.t labels], reads);
    saveas(h.hgAxes, sprintf('pictures/S%s%d_particle_%d.jpg', strSamples, job_id, i));
end

%%
nClasses = 3;
priors.R = generate_sticky_prior(0.0001, 4, 630000);
priors.NT = generate_sticky_prior(0.001, 4, 630);
priors.AA = generate_sticky_prior(1/21, 21, 210);
priors.decay = [0 1e-3];  % mean and std
priors.pi = 100*ones(1,64);
priors.nClasses = 3;

%                      generate_data(  str,  nReads, L, maxCells, maxTime, nClasses, br,  dr)       
[reads labels truth] = generate_data('synth', 1000, 33, 4200,     10,      nClasses, 0.8, 0.5, priors);


[a stats stats2] = infer_tree(reads, 1, priors);





end

