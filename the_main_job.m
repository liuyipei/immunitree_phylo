% Experiment running
run_as_a_job = true;
synthetic = false;
%% path
cd ~/JVL/src/phylo/
addpath('../phylo');
addpath('../DPtrees/');
addpath('../VDJ/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));
global codon2aa codon2nt nt2codon; [codon2aa codon2nt nt2codon]= get_maps();

RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',sum(100*clock)));

dir = '~/scratch/data/Fire_FL'; 

if ~exist('run_as_a_job', 'var')
%% starting point:  Load/generate reads
    clear; synthetic_DB = {'synthetic_test_1_MAY_11_2011', 'synthetic_02_23_11'};
    synthetic = true;  samples = synthetic_DB{1};
    %%
    clear; synthetic = true;  samples = [];
    %%
    clear; synthetic = false;  samples = [68 69];
    %%
    clear; synthetic = false;  samples = 'just_in_result_72_73_85388';
end
%%
synthetic, samples
close all;
job_id = ceil(rand*100000);
priors = [];
filename = [];
is_codon = false;

if synthetic
    if ~isempty(samples)
        assert(ischar(samples));
        load(samples);
    else        
        priors_gen = struct('is_codon', is_codon, 'is_decay', false, 'nClasses',3);
        %                      generate_data(  str,  nReads, L, maxCells, maxTime, br,  dr)       
        [reads labels truth] = generate_data('synth', 1000, 100, 4200,     10,      0.8, 0.5, priors_gen);
        filename = sprintf('synthetic_%d', job_id);
        priors.pi = truth.sequences(1,:);
    end
else % real   
    if ischar(samples)
        load(samples);
        samples = truth.samples;
    end
    [reads labels truth] = generate_data('real', samples, is_codon); % false for nt mode
    truth.labels = labels; labels = [];
    priors.pi = truth.pi;
    strSamples = sprintf('%d_', samples);
    filename = sprintf('S%s%d', strSamples, job_id);
end
    

%% visualize ground truth (erasing dead subtrees)
if isfield(truth, 'tree')
    a = convert_phylo_tree_to_mutation_tree(truth);
    if is_codon
        visualize_tree(a.tree, codons2seqs(a.sequences), [a.t labels], reads);
    else
        visualize_tree(a.tree, a.sequences, [a.t labels], reads);
    end
end


%% init MCMC (all we need is 'reads').

priors.is_codon = is_codon;
priors.is_decay = false;
priors.nClasses = 1;
%priors.R = generate_sticky_prior(1e-7, 4, 6300000);

nIter = 1000;
[best cur stats stats2] = infer_tree(reads, nIter, priors, filename, labels);


if exist('run_as_a_job', 'var')
    if length(samples) == 1
        fid = fopen(sprintf('%s/jobs/S%d.done', dir, samples), 'w');
        fprintf(fid, sprintf('S%s%d\n', strSamples, job_id));
        fclose(fid);
    end
    exit
end



%%

if 0

clear 
is_codon = false;
samples = [70 71];
job_id = 75239;


%%
strSamples = sprintf('%d_', samples);
filename = sprintf('S%s%d', strSamples, job_id);

[reads labels truth] = generate_data('real', samples, is_codon); % false for nt mode
load(sprintf('S%s%d.mat', strSamples, job_id));
   
% visualize best reconstruction (erasing dead subtrees)
a = best;
a = convert_phylo_tree_to_mutation_tree(a);
a = collapse_edges(a);
if is_codon
    h = visualize_tree(a.tree, codons2seqs(a.sequences), [a.t labels], reads);
else
    h = visualize_tree(a.tree, a.sequences, [a.t labels], reads);
end

%
figure(4);
if isfield(truth, 't')                                  
    report_scores(truth, a, labels);
else
    report_scores(reads, a);
end
%
report_traces(best, stats2, stats);


%%  save figures to file
figure_positions = [0           0           0           0;
        1084         645         721         165;
        607         544        1021         953
        1815         645         681         165;
        1404         894        1088         603];

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
    

end


