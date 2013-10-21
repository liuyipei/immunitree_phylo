% input:
% T is Mx1 array.  T(i) is the parent of node i.
% singles is MxBxL, where L is the number of sites (length of each seq).
%    At the end this will store the marginal probabilities for each node of
%    the tree.  Now it stores some prior information on each node.
% pairwise: a BxBxL array with pairwise 
%    potentials for each edge in the tree. pairwise(i,j,l) is the probability
%    a parent with base j generated a child with base i, on location l.
%    exception:  If pairwise is a BxB matrix, then convert to BxBxL by:
%        Q = Q';  %reverse parent_child to child_parent potential
%        Q = Q(:,:,ones(1,L));
function [sequences prob] = tree_sample_for_phylo_unlog(T, singles, pairwise, seq_root_parent, sequences, annealing_prob, seq_root_parent_weight)

global greedy;
if isempty(greedy), greedy = 0; end
if nargin == 0, sequences = test_tree_sample(); return; end
% the weight of the prior sequence should be adjusted downwards, to less
% than 1. (e.g. 0.5 ~ given weight as half of a data point)
if ~exist('seq_root_parent_weight', 'var'), seq_root_parent_weight = 0.5; end 
if ~exist('annealing_prob', 'var'), annealing_prob = 0; end
annealed = rand(1) < annealing_prob;

B = size(pairwise, 1);
   
M = length(T);

if isscalar(singles) % no evidence, sampling from prior. singles is the length of each seq
    evidence = false;
    L = singles;
else
    evidence = true;
    L = size(singles, 3);
end

pairwise = reshape(pairwise,B,B*L); 

% Phase one - compute unnormalized conditionals:
% go over tree in reverse topological ordering
DG = sparse(1+T, 2:M+1, true, M+1, M+1);
if M == 1, 
    order = 1;
else
    order = graphtopoorder(DG(2:end, 2:end));
end
rev_order = fliplr(order);

ix = 1:L;
ix = ix(ones(1,B), :);
ix = ix(:);

pairwise_1bp = pairwise(1:B, 1:B); % 5x5 transition matrix for yi's optimization
yi_debug = 0; % check correctness of yi's version
if evidence
    for k = rev_order(1:end-1)
        % each nodes updates the singleton potential of its parent
        % hence, once we get to a node, it was updated by all its children
        % ==> its singleton = its conditional
        % multipile my own conditional with pairwise-potential
        % child_parent(i,j,l) = pairwise_k(i,j,l) + singleton_k(j, l)
        
        if yi_debug
            % yi's optimized implementation -- check it
            yi_message=reshape(pairwise_1bp * reshape(singles(k,:, :), B, L), 1, B*L);
            %joni's sensible implementation
            child_parent = reshape(repmat(reshape(singles(k,:, :), B, L), B, 1), B, B*L) .* pairwise;

            % TODO:  The above line is expensive.  But I can't make it any faster.
            %child_parent_ = squeeze(singles(k,:, ix)) .* pairwise;
            %assert(max(abs(child_parent(:)-child_parent_(:)))<1e-10);

            % send message to parent - update singleton of parent
            % child_parent is 4x4L.  Sum over child to infer parent.
            message = sum(child_parent);  % message is 1x4L
            assert(max(abs(yi_message - message))<12*eps);
        else
            % just do yi's optimized implementation
            message=reshape(pairwise_1bp * reshape(singles(k,:, :), B, L), 1, B*L);
        end
        
        %message = message./sum(message);  % normalize for numerical stability
        
        singles(T(k), :, :) = singles(T(k), :, :) .* reshape(message, 1, B, L); % 1x4xL            
        partition = sum(singles(T(k), :, :));
        singles(T(k), :, :) = singles(T(k), :, :) ./ partition(1,ones(1,B), :);
    end    
end

if ~exist('sequences', 'var'), sequences = greedy*ones(M,L); end
if isempty(sequences), sequences = greedy*ones(M,L); end

sequences = double(sequences);
assert(isequal(size(sequences), [M L]));

prob = zeros(M,L);
cpd = 1;
muls = (0:B:(B*L-B));

% handle root
k = order(1);
assert(T(k) == 0);
if evidence, cpd = reshape(singles(k,:, :), B, L); end

if nargin<4 || isempty(seq_root_parent)
    % sample root bases from uniform prior
    [sequences(k,:) prob(k,:)] = many_mnrnd(cpd .* ones(B,L), sequences(k,:));
else
    if size(seq_root_parent,1) == B % discrete distribution (in log-space) is given over root
        if size(seq_root_parent,2) ~= L
            seq_root_parent = seq_root_parent(:,ones(1,L));
        end
        % sample root from given prior (plus likelihood)
        [sequences(k,:) prob(k,:)] = many_mnrnd(cpd .* seq_root_parent, sequences(k,:));
    else
        % root seq completely given.
        
        % last swapped August 7, 2012
        % root_pairwise = [generate_sticky_prior(0.01, B, 1) (1/B)*ones(B,1)]; % from germline to clone-root
        
        % use average of pairwise mutation model and uniform model        
        root_pairwise = ((1-seq_root_parent_weight) * pairwise_1bp + seq_root_parent_weight * ones(B,B)/B) / 2;
        if annealed            
            [prob(k,:) sequences(k,:)] = max(cpd .* root_pairwise(:, seq_root_parent), [], 1);
        else % sample
            [sequences(k,:) prob(k,:)] = many_mnrnd(cpd .* root_pairwise(:, seq_root_parent), sequences(k,:));
        end
    end
end        

given = (sequences(order,1) > 0);
all_given = (sum(~given) == 0);

muls_BB= 0:(B*B):(B*B*(L-1));
if annealed
    if ~evidence
        if ~all_given            
            for k = order(~given)
                parent_ix = sequences(T(k), :) + muls;
                [prob(k,:) sequences(k,:)] = max(pairwise(:,parent_ix), [], 1);
            end
        end
    else
        for k = order(2:end)
            % this loop is the most expensive in the program.
            cpd = reshape(singles(k,:, :), B, L);
            parent_ix = sequences(T(k), :) + muls;
            [prob(k,:) sequences(k,:)] = max(cpd .* pairwise(:, parent_ix), [], 1);
        end
    end
else
    if ~evidence
        if ~all_given
            assert(greedy == 0);
            cdf = cumsum(pairwise,1);
            assert( abs(cdf(end,1)-1)<1e-5);
            rnd = rand(size(sequences));    
            for k = order(~given)
                parent_ix = sequences(T(k), :) + muls;
                sequences(k,:) = sum(cdf(:,parent_ix) < rnd(k*ones(1,B), :), 1) + 1;
            end
        end
        if nargout > 1
            % calculate mutation probabilities faster given all the sequences
            non_root = [1:(order(1)-1) (order(1)+1):M];
            parent = sequences(T(non_root,1),:);
            child = sequences(non_root,:);
            ix = mysub2ind(B, child(:), parent(:), 1);
            ix = reshape(ix,M-1,L);            
            ix = ix + muls_BB(ones(M-1,1),:);        
            prob(non_root,:) = pairwise(ix);
        end        
    else
        for k = order(2:end)
            % this loop is the most expensive in the program.
            cpd = reshape(singles(k,:, :), B, L);
            parent_ix = sequences(T(k), :) + muls;
            [sequences(k,:) prob(k,:)] = many_mnrnd(cpd .* pairwise(:, parent_ix), sequences(k,:), evidence);        
        end          
    end
end
sequences = int16(sequences);

end


% TODO:  unit-test

function sequences = test_tree_sample()
%%
epsilon = 0.01;
sigma = 0.05;

T     = [0 1 1 2 2 3 3]';               %     1
                                        %  1    2
truth = [1 1 2 1 3 2 4]';               % 1 3  2 4
miss  = [2 2 3 2 4 3 2]';

counts= [4 3 4 3 4 3 3]';

M = length(T);
singles = zeros(M, 4);
for k=1:length(T)
    singles(k,truth(k)) = (counts(k)-1)*(log(1-epsilon)-log(epsilon));
    singles(k,miss(k))  = log(1-epsilon)-log(epsilon);
end
exp(singles)
pairwise = (1-4/3*sigma)*eye(4) + sigma/3*ones(4)
[sequences prob] = tree_sample_for_phylo_unlog(T, exp(singles), pairwise)
truth
assert(isequal(truth, sequences));
fprintf('Test Passed\n');
%%
end

function res = log_sum_exp(v)
    v_max = max(v);
    v = v-v_max(ones(1,size(v,1)), :);
    res = log(sum(exp(v))) + v_max;
end


function sequences = test_tree_sample2()
%%
epsilon = 0.01;
sigma = 0.05;

T     = [0 1 2];
truth = [1 1 1];
counts= [100 0 1];

M = length(T);
singles = zeros(M, 4);
for k=1:length(T)
    singles(k,truth(k)) = counts(k)*(log(1-epsilon)-log(epsilon));
end
exp(singles)
pairwise = (1-sigma-sigma/3)*eye(4) + sigma/3*ones(4)
sequences = tree_sample(T, singles, log(pairwise))'
truth
%%
end




