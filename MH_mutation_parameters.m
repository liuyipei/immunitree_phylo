% function [AA NT pi codes total_mutations ll] = MH_mutation_parameters(tree, codons, AA, NT, rate_class, priors)
% updates the mutation potentials AA,NT,and pi using an MH move.
% returns the total number of mutations in the current configuration, as
% well as the log likelihood over the value of the potential and the
% sequences.
function [mut_model codes total_mutations ll] = MH_mutation_parameters(tree, codons, t, mut_model, priors)
%   it is no longer true that we can gather the class statistics from all
%   the sites at once, because the decay factor is going to make the
%   transition matrices different.  


assert(isequal(find(tree(:,1) == 0),1));
is_codon = priors.is_codon;

% work only on the subtree tied to evidence
evidence = nodes_with_evidence(tree,t);
tree(~evidence,1) = -1;
[tree, codons] = clean_tree(tree,codons);

total_mutations = sum(sum(codons(2:end,:) ~= codons(tree(2:end,1),:)));

codes = zeros(3,2);    

% compute pi (prior over root) (for synthetic data only)
if size(priors.pi,2) == size(codons,2)
    mut_model.pi = priors.pi;
else
    mut_model.pi = drchrnd(histc(codons(1,:), 1:size(mut_model.pi,2)) + priors.pi,1);
end

% if called only to compute likelihood then:
if nargout == 2 % [ll total_mutations] = ...
    [l_data l_prior] = get_seq_likelihood(tree, codons, mut_model, priors);
    codes = total_mutations;
    mut_model = [sum(l_data) sum(l_prior)];
    return;
end

if is_codon
    % sample AA + compute likelihood of AA prior
    [mut_model.AA codes(1,:)] = MH_AA(tree, codons, mut_model, priors);
end

% sample NT + compute likelihood of sequences and NT prior
[mut_model.NT codes(2,:) l_data l_prior] = MH_NT(tree, codons, mut_model, priors);

% sample decay parameter - from prior
if sum(priors.decay) > 0
    [mut_model.decay l_data codes(3,:)] = MH_decay(tree, codons, mut_model, priors, l_data);
end

% combine likelihoods to the returned variable ll
ll = [sum(l_data) sum(l_prior)];

end

function [l_data l_prior] = get_seq_likelihood_AA(tree, codons, mut_model, priors)
    global codon2aa
    % compute likelihood of sequences    
    B = size(mut_model.AA,1);
    mut_model.AA = reshape(mut_model.AA, B, B, []);
    nClasses = size(mut_model.AA,3);
    Q = get_Qs(mut_model);
    [~,prob ] = tree_sample_for_phylo_unlog(tree(:,1), size(codons,2), Q, mut_model.pi', codons);
    prob = log(prob(2:end,:));

    
    % create l_data, l_prior:  matrices of size 21 x nClasses
    % l_data(aa,l) = total probability for mutations from aa to
    %                something else, in class l.
    l_data = zeros(B, nClasses);
    l_prior = zeros(B, nClasses);
    parents = codon2aa(codons(tree(2:end,1), :));    
%    rate_class = rate_class(ones(1,size(parents,1)),:);
    for l = 1:nClasses
        iz = (mut_model.rate_class == l);
        prob_ = prob(:,iz);
        parents_ = parents(:,iz);
        for aa = 1:B            
            l_data(aa,l) = sum(prob_(parents_ == aa));
        end
        [~, l_prior(:,l)] = my_gampdf(mut_model.AA(:,:,l), priors.AA);
    end
    l_data = l_data(:)';
    l_prior = l_prior(:)';
end


% computes likelihood of sequences and likelihood of NT and AA
function [l_data l_prior] = get_seq_likelihood(tree, codons, mut_model, priors)
    B = 5;
    % compute likelihood of sequences    
    mut_model.NT = reshape(mut_model.NT, B, B, []);
    Q = get_Qs(mut_model);        
    nClasses = size(mut_model.NT,3);
    
    [~,prob] = tree_sample_for_phylo_unlog(tree(:,1), size(codons,2), Q, mut_model.pi', codons);
    prob = sum(log(prob(2:end,:)),1); % prob is 1xL
    % TODO: if priors.is_codon == false, we might want to break prob
    % according to parents, similar to get_seq_likelihood_AA

    
    % create l_data, l_prior:  matrices of size nClasses x 21
    % l_data(l,aa) = total probability for mutations from aa to
    %                something else, in class l.
    l_data = zeros(1, nClasses);
    l_prior = zeros(1, nClasses);
    for l = 1:nClasses
        iz = (mut_model.rate_class == l);
        l_data(l) = sum(prob(iz));

        l_prior(l) = 0;
        if nargout > 1
            l_prior(l) = l_prior(l) + my_gampdf(mut_model.NT(:,:,l), priors.NT);
            if priors.is_codon
                l_prior(l) = l_prior(l) + my_gampdf(mut_model.AA(:,:,l), priors.AA);
            end
        end
    end
end


function [AA codes l_prior] = MH_AA(tree, codons, mut_model, priors)
    codes = [0 0];
    
    % l_data = the likelihood of the sequences determined by each column of
    % AA.  l_prior = the prior likelihood over each column of AA.
    [l_data l_prior] = get_seq_likelihood_AA(tree, codons, mut_model, priors);
    
    AA = reshape(mut_model.AA, size(mut_model.AA,1), []);
        
    % take one coordinate in every column, and chage its weight    
    for i=1:size(AA,1)
        r = randn;
        AA_ = AA;
        AA_(i,:) = AA_(i,:)*exp(r);
        
        % calculate likelihood
        mut_model.AA = AA_;
        [l_data_ l_prior_] = get_seq_likelihood_AA(tree, codons, mut_model, priors);
        
        % transitions will cancel accept Jacobian term = exp(r)
        acceptance = exp(l_data_-l_data + l_prior_-l_prior)*exp(r); % size = [1, 21*nClasses]
        accepted = rand(1,size(AA,2)) < acceptance; 
        AA(i,accepted) = AA_(i,accepted);
        l_data(accepted) = l_data_(accepted);
        l_prior(accepted) = l_prior_(accepted);
        codes = codes+[sum(accepted) size(AA,2)];              
    end
    
    AA = reshape(AA, size(AA,1), size(AA,1), []);
        
end

function [NT codes l_data l_prior] = MH_NT(tree, codons, mut_model, priors)
    codes = [0 0];
    B = size(mut_model.NT,1);
    
    % l_data = the likelihood of the sequences determined by each column of
    % AA.  l_prior = the prior likelihood over each column of AA.
    [l_data l_prior] = get_seq_likelihood(tree, codons, mut_model, priors);
    
    % TODO:  In case priors.is_codon == false, this could be done more
    %        efficiently, by sampling one entry from every column of NT.
    %        All you need to change here is the reshape to:
    %        NT = reshape(mut_model.NT, B, []); % NT is 4 x 4*nClasses
    %        and then change the function get_seq_likelihood accordingly
    NT = reshape(mut_model.NT, B*B, []); % NT is 16 x nClasses
    for i=1:size(NT,1)
        r = randn;
        NT_ = NT;
        NT_(i,:) = NT(i,:)*exp(r);

        % calculate likelihood
        mut_model.NT = NT_;
        [l_data_ l_prior_] = get_seq_likelihood(tree, codons, mut_model, priors);
        
        % transitions will cancel accept Jacobian term = exp(r)
        acceptance = exp(l_data_-l_data + l_prior_-l_prior)*exp(r); % 1 x nClasses
        accepted = rand(1,size(NT,2))<acceptance;
        NT(i,accepted) = NT_(i,accepted);
        l_data(accepted) = l_data_(accepted);
        l_prior(accepted) = l_prior_(accepted);
        codes = codes+[sum(accepted) size(NT,2)];              
    end
            
    NT = reshape(NT, B, B, []);
        
end


function [ll ll_decomposed] = my_gampdf(phi, prior)
    if isstruct(prior)        
        epsilon = prior.eps; pseudo = prior.pseudo; B = prior.dim;
        a = pseudo*(epsilon + (1-B*epsilon)*eye(B));
    else
        a = prior;
    end
    ll_decomposed = sum((a-1).*log(phi)) -sum(phi) -sum(gammaln(a));
    ll = sum(ll_decomposed);
end



function [decay l_data codes] = MH_decay(tree, codons, mut_model, priors, l_data_prev)
    codes = [0 1];
    decay = mut_model.decay;
    mut_model.decay = mut_model.decay + randn*(0.5*priors.decay(2));
    Q = get_Qs(mut_model);        
    [~,prob] = tree_sample_for_phylo_unlog(tree(:,1), size(codons,2), Q, mut_model.pi', codons);
    l_data = sum(sum(log(prob(2:end,:)))); 
%    l_data = get_seq_likelihood(tree, codons, rate_class, pi, AA, NT, decay_);

    log_decay_prior = @(x, prior) -(x-prior(1) )^2 / (2*prior(2)^2);
    prior_diff = log_decay_prior(mut_model.decay, priors.decay)-log_decay_prior(decay, priors.decay);
    if rand < exp(sum(l_data)-sum(l_data_prev) + prior_diff)
        decay = mut_model.decay;
        codes(1) = 1;
    end    
end