function [l_data l_prior] = get_seq_likelihood(SS, AA, NT, AA_prior, NT_prior, decomposed)

    % computation
    Q = make_codon_matrix(AA, NT);

    if ~decomposed  
        l_prior = my_gampdf(AA, AA_prior) + my_gampdf(NT, NT_prior);
        l_data = SS(:)'*log(Q(:));
    else % decomposed to AA
        % working on AA 
        global codon2aa    
        [~, l_prior] = my_gampdf(AA, AA_prior);
        l_data = sum(SS.*log(Q));
        l_data = l_data*sparse(1:64, double(codon2aa), true);
    end
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

