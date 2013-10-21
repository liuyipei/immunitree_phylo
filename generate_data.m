function [reads labels truth] = generate_data(str, nReads, L, maxCells, maxTime, birth_rate, death_rate, priors)

if strcmp(str, 'real')    
    if ~exist('nReads', 'var'), assert(false); end;
    samples = nReads; 
    [Z labels germline] = parse_dna(samples);
    Z = int16(Z);
    [I, ~] = find(Z == 5);
    Z(I, :) = [];
    reads = int16(Z);
    labels(I) = [];
    truth.pi = germline.seq;
    if L == true % codon mode
        trim_left = mod(3-germline.frame_shift, 3);
        trim_right = mod(size(reads,2)-trim_left, 3);
        reads = reads(:,(1+trim_left):(end-trim_right));
        truth.pi = seqs2codons(germline.seq((1+trim_left):(end-trim_right)));
    end
    truth.desc = ['Real cancer data. samples: ' sprintf('%d ', samples)];
    truth.samples = samples;
    labels = labels';
    
elseif strcmp(str, 'synth')
    
    rates = struct('birth', birth_rate, 'death', death_rate);
    priors = init_priors(priors);
    
    % generate the tree
    live_cells = [];
    while isempty(live_cells)
        [tree, sequences, F, mut_model] = ...
            generate_sequences(L, maxCells, maxTime, rates, priors);
        live_cells = find(tree(:,3) == F);
        fprintf('Total %d live cells from %d generated cells. Time stopped at F=%.2f.\n', ...
            length(live_cells), size(tree,1), F);
    end

    % Now generate reads from live cells
    t = live_cells(ceil(length(live_cells)*rand(1,nReads)));

    % read model
    R = generate_stochastic_matrix(priors.R, 1, false);
    reads = sequences(t,:);  % assume no noise for now.
    if priors.is_codon
        reads = codons2seqs(reads); % the nucleotide version of sequences
    end

    for i=1:length(t)        
        reads(i,:) = tree_sample_for_phylo(0, size(reads,2), eye(4), R(reads(i,:),:)');
    end
    
    truth.desc = 'synthetically generated data';
    truth.samples = '';
    truth.tree = tree;
    truth.sequences = sequences; 
    truth.t = t; 
    truth.R = R; 
    truth.rates = rates;
    truth.F = F;
    truth.mut_model = mut_model;

    labels = convert_phylo_tree_to_mutation_tree(truth);
    labels = labels.t;
else
    assert(false);
end

end % function

function [tree, sequences, F, mut_model] = ...
    generate_sequences(L, S, F, rates, priors)

nClasses = priors.nClasses;

if priors.is_codon
    L = floor(L/3);
    mut_model.AA = generate_stochastic_matrix(priors.AA, nClasses);
end

mut_model.decay = priors.decay(1)+randn*priors.decay(2);
mut_model.NT = generate_stochastic_matrix(priors.NT, nClasses);
mut_model.pi = drchrnd(priors.pi, 1);
mut_model.nClasses = nClasses;
%mut_model.B = 64;  % no one is using this field
%mut_model.L = L;    % is anyone using this field?

% decide for each location, if it is conserved or not
logphi = log(drchrnd(ones(1, nClasses),1));
mut_model.rate_class = many_lmnrnd(logphi(ones(L,1), :)', 0);

% mutation model

% TODO: replace this code with 
Qs = get_Qs(mut_model);
[tree, sequences, F] = generate_sequences_from_prior(S, F, rates, Qs, mut_model.pi');


end



function mut_rates = get_mutation_rates()
  assert(false);  % not used!
% How to compute 3 classes of mutation rates, in the case the rate is
% distributed according to gamma(1,1):
    mut_rates = zeros(1,3);
    dx = 0.01; 
    x = dx:dx:gaminv(1/3,1,1);  
    mut_rates(1) = 3*x*dx*gampdf(x, 1, 1)';
    x = gaminv(1/3,1,1):dx:gaminv(2/3,1,1);  
    mut_rates(2) = 3*x*dx*gampdf(x, 1, 1)';
    x = gaminv(2/3,1,1):dx:10;  
    mut_rates(3) = 3*x*dx*gampdf(x, 1, 1)';
end

