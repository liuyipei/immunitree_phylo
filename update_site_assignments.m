%function [rate_class, ll] = update_site_assignments(rate_class, AA, NT, decay, pi, tree, sequences)
function [rate_class, ll] = update_site_assignments(mut_model, tree, sequences)

    global greedy
    rate_class = mut_model.rate_class;
    
    % updates mutation classes
    % assuming frequencies were drawn from a Dirichlet
    K = mut_model.nClasses; 
    L = length(mut_model.rate_class);
    log_prob = zeros(K,L);
    for k=1:K
%        Qs = get_Qs(Q, k*ones(1,L));
        mut_model.rate_class = k*ones(1,L);
        Qs = get_Qs(mut_model);
        [~,prob] = tree_sample_for_phylo_unlog(tree(:,1), L, Qs, mut_model.pi', sequences);        
        log_prob(k,:) = sum(log(prob));
    end
    
    pseudo = 1;
    N = hist(rate_class, 1:K)+pseudo;    
    for l=randperm(L)
        N(rate_class(l)) = N(rate_class(l))-1;
        rate_class(l) = lmnrnd(log_prob(:,l)'+log(N), greedy);
        N(rate_class(l)) = N(rate_class(l))+1;
    end
    assert(isequal(N, hist(rate_class, 1:K)+pseudo)); 
    
    ll = sum(gammaln(N)-gammaln(pseudo)) + gammaln(K*pseudo) - gammaln(sum(N));

%     muls = 0:K:(L*(K-1));
%     ll_seqs = log_prob(muls+rate_class);
    

end
