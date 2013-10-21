%function [tree sequences t F rates mut_model R] = ...
%    initialize_chain(reads, leaf_level, priors, F, init_)
% set the initial state of the Markov Chain
% Inputs: 
%  reads: the data.
%  leaf_level: code for different initializations methods for the tree.
%              surrently supports only leaf_level = 0.
%  priors: priors on probabilistic model.  Most are documented under
%          'init_priors'.  Note that priors.pi can be set either to the 
%          Dirichlet weights for the prior over the root sequence, OR
%          it can be set as the root sequence itself (in that case the root
%          sequence is not sampled, except from the locations priors.pi has 
%          gaps).
%  F: snapshot time
%  init_: a state given by the user with some of the variable values given.
%         the function fills in only the missing values of init_.
%
% Outputs are the state variables:
%  F: the snapshot time. Cells alive at this time are considered "alive".
%  tree: each row represents a single cell in the tree, and states its
%         parent node, birth time, and death time. If the cell is alive,
%         its death_time = F.
%  sequences: matrix of all cell sequences, ordered by their cell index.
%  t: assignments of reads to cells in the tree.  
%  rates: a structure with birth and death rates for tree.
%  mut_model: a structure with all the parameters governing the mutatio
%             model.
%  R: (used to be) a 5x4 transition matrix for the read noise model. 
%      the last row is all ones, it's there to support reads with gaps.
%     (as of Yi) now a 5x5 proper transition matrix, with gap being a real
%      character
function [tree sequences t F rates mut_model R] = ...
    initialize_chain(reads, leaf_level, priors, F, init_)

nClasses = priors.nClasses;
is_codon = priors.is_codon;
L = size(reads,2);
%init.rates.birth = 1;
init.rates.birth = 0.5; % changed from 1 to 0.5 in August 2008
init.rates.death = 0.5;
init.mut_model.decay = 0;
init.mut_model.nClasses = nClasses;        
init.mut_model.NT = generate_stochastic_matrix(priors.NT, nClasses);
if is_codon
    init.mut_model.AA = generate_stochastic_matrix(priors.AA, nClasses);
end
if is_codon
    assert(mod(L, 3)== 0);
    init.mut_model.rate_class = ceil(nClasses*rand(1, L/3));
else
    init.mut_model.rate_class = ceil(nClasses*rand(1, L));
end
init.R = [generate_stochastic_matrix(priors.R, 1); ones(1,5)];
    

if leaf_level==0
    % create a tree where every set of identical reads belong to a node
    % that is a direct child of the root

    % Fill-in the gaps in the reads when you generate the sequences
    % because we don't allow gaps ('N') in the sequences
    sequences = reads;
    consensus = histc(reads, 1:5, 1);
    [~,consensus] = max(consensus(1:4,:));
    consensus = consensus(ones(1,size(sequences,1)),:);
    sequences(sequences == 5) = consensus(sequences == 5);        
    
    % For each set of unique reads, set one tree cell with that sequence
    % and assign the set's reads to that cell
    [sequences, ~, t] = unique(sequences, 'rows'); 

    % create the tree where all those cells are direct children of the root
    tree = zeros(size(sequences,1), 3);
    tree(:,3) = F;
    tree(:,2) = F/2 * rand(size(tree,1),1);
    tree(:,1) = 1;
    tree = [0 0 F/2 ; tree];
    root_seq = median(double(sequences),1);
    sequences = [root_seq; sequences]; t = t+1;


else % find minimum spanning tree.  
    assert(false); % Functionality not supported
end

if is_codon % switch sequences in tree to codons.
    sequences = seqs2codons(sequences);
end
sequences = int16(sequences);

% two cases for priors.pi:
if size(priors.pi,2) == size(sequences,2)
    % priors.pi is the root sequence itself.
    init.mut_model.pi = priors.pi;
else
    % priors.pi is the dirichlet prior weights over root sequence letters.
    % initialize the root by sampling from this prior.
    init.mut_model.pi = drchrnd(priors.pi, 1);
end

init.F = F;
init.t = t;
init.tree = tree;
init.sequences = sequences;

% override elements of init state set by the user:
% First, override elements of mutation model set by the user
if isfield(init_, 'mut_model')
    init_.mut_model = struct_override(init.mut_model, init_.mut_model);
end
init = struct_override(init, init_); % then the rest of the state

[tree sequences t F rates mut_model R] = deal(init.tree, init.sequences, ...
    init.t, init.F, init.rates, init.mut_model, init.R);

end


% % Unsupported functionality for creating minimum spanning tree over the set
% % of reads
%     N = size(reads,1);
%     root_seq = median(reads,1);
%     sequences = [root_seq; reads];
%     map = 'ACGT';
%     fprintf('gathering distances...');
%     D = pdist(double(sequences), 'hamming')*size(sequences,2);
%     %D = seqpdist(map(sequences),'method','jukes-cantor','Alphabet', 'NT');
%     fprintf('Done.\n');
% 
%     % choose a root r
%     fprintf('computing MST...');
%     [~, pred] = graphminspantree(sparse(squareform(D)), 1);
%     fprintf('Done.\n');
%     pred = pred';
%     assert(pred(1) == 0);
%     t = 2:N+1;
%     for i=1:leaf_level%
%         has_children = zeros(1,length(pred));
%         has_children(pred(2:end,1)) = 1;
%         leaves = find(~has_children);
%         reads_in_leaves = find(ismember(t,leaves));
%         t(reads_in_leaves) = pred(t(reads_in_leaves));
% %        t(leaves-1) = pred(leaves);
%         pred(leaves) = -1;
%         [pred, sequences, t] = clean_tree(pred, sequences, t);
%     end
% 
%     M = length(pred);
%     tree = [pred zeros(M,2)];  
%     assert(tree(1,1) == 0);
%     DG = sparse(1+tree(:,1), 2:M+1, true, M+1, M+1);
%     %DG = get_graph_from_tree(T, true, true); % should be the same as above
% 
%     order = graphtopoorder(DG(2:end, 2:end));
%     for k=order(2:end)  % needs to be topological order
%         tree(k,2) = tree( tree(k,1),2) + exprnd(1);
%     end
%     F_ = max(tree(:,2))+1;
%     tree(:,2) = tree(:,2)*F/F_;
%     tree(:,3) = F;
%     ix = find(tree(:,1) == 1);
    


