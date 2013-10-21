% input:
% T is Mx1 array.  T(i) is the parent of node i.
% singles is Mx4xL, where L is the number of sites (length of each seq).
%    At the end this will store the marginal probabilities for each node of
%    the tree.  Now it stores some prior information on each node.
% pairwise: a BxBxL array with pairwise 
%    potentials for each edge in the tree. pairwise(i,j,l) is the probability
%    a parent with base j generated a child with base i, on location l.
%    exception:  If pairwise is a BxB matrix, then convert to BxBxL by:
%        Q = Q';  %reverse parent_child to child_parent potential
%        Q = Q(:,:,ones(1,L));
function [sequences prob] = tree_sample_for_phylo(T, singles, pairwise, seq_root_parent, sequences)

global greedy;
if isempty(greedy), greedy = 0; end
if nargin == 0, sequences = test_tree_sample(); return; end

B = size(pairwise, 1);
   
M = length(T);

evidence = true;
if isscalar(singles) % no evidence, sampling from prior. singles is the length of each seq
    evidence = false;
    singles = zeros(M,B,singles);  % size: [M,B,L]
end
L = size(singles, 3);

Q = pairwise;  % Q is the pairwise potential in log scale.
if length(size(Q)) == 2 && size(Q,1) == size(Q,2)
    Q = Q';  %reverse parent_child to child_parent potential
    Q = Q(:,:,ones(1,L));
end    
Q = reshape(Q,B,B*L); 
pairwise = Q;    
clear Q;

% Phase one - compute unnormalized conditionals:
% go over tree in reverse topological ordering
DG = sparse(1+T, 2:M+1, true, M+1, M+1);
if M == 1, 
    order = 1;
else
    order = graphtopoorder(DG(2:end, 2:end));
end
rev_order = fliplr(order);

if evidence
    for k = rev_order
        % each nodes updates the singleton potential of its parent
        % hence, once we get to a node, it was updated by all its children
        % ==> its singleton = its conditional
        if T(k) ~= 0 % if there is a parent
                % multipile my own conditional with pairwise-potential
                % child_parent(i,j,l) = pairwise_k(i,j,l) + singleton_k(j, l)
                child_parent = reshape(repmat(reshape(singles(k,:, :), B, L), B, 1), B, B*L) + pairwise;

                % send message to parent - update singleton of parent
                % child_parent is 4x4L.  Sum over child to infer parent.
                message = log_sum_exp(child_parent);  % message is 1x4L
                singles(T(k), :, :) = singles(T(k), :, :) + reshape(message, 1, B, L); % 1x4xL            
        end
    end
end

if ~exist('sequences', 'var'), sequences = greedy*ones(M,L); end
assert(isequal(size(sequences), [M L]));

prob = zeros(M,L);

for k = order
    log_cpd = reshape(singles(k,:, :), B, L);
    if T(k) == 0 % if there is no parent
        if nargin<4 || isempty(seq_root_parent)
            % sample root bases from uniform prior
            [sequences(k,:) prob(k,:)] = many_lmnrnd(log_cpd, sequences(k,:));        
        else
            if size(seq_root_parent,1) == B % discrete distribution (in log-space) is given over root
                if size(seq_root_parent,2) ~= L
                    seq_root_parent = seq_root_parent(:,ones(1,L));
                end
                % sample root from given prior (plus likelihood)
                [sequences(k,:) prob(k,:)] = many_lmnrnd(log_cpd + log(seq_root_parent), sequences(k,:));
            else % root seq completely given.  set root as given. (sample uniformly missing letters).
                root_pairwise = log([eye(B) (1/B)*ones(B,1)]);
                [sequences(k,:) prob(k,:)] = many_lmnrnd(log_cpd + root_pairwise(:, seq_root_parent), sequences(k,:));
            end
        end        
    else      
        % sample node from given parent 
        parent_ix = double(sequences(T(k), :)) + (0:B:(B*L-B));
        [sequences(k,:) prob(k,:)] = many_lmnrnd(log_cpd + pairwise(:, parent_ix), sequences(k,:));
    end
end

end


% TODO:  unit-test

function sequences = test_tree_sample()
%%
epsilon = 0.01;
sigma = 0.05;

T     = [0 1 1 2 2 3 3]';

truth = [1 1 2 1 3 2 4]';
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
[sequences prob]= tree_sample_for_phylo(T, singles, log(pairwise))
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




