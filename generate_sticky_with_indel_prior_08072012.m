function prior = generate_sticky_with_indel_prior_08072012(epsilon, B, pseudo, perc_mut_indel)
    %prior = pseudo * (epsilon + (1-B*epsilon)*eye(B));
    prob_point_mut = (1-perc_mut_indel)*epsilon/(B-2);
    prob_indel = epsilon * perc_mut_indel;
    
    prior = ones(B) - eye(B);
    
    prior(1:(B-1), 1:(B-1)) = prior(1:(B-1), 1:(B-1)) .* prob_point_mut;
    
    prior(B, (1:B-1)) = prob_indel;
    prior(1:(B-1), B) = prob_indel;
    
    for i = 1:B
        prior(i,i) = 0;
        prior(i,i) = 1 - sum(prior(i,:));
    end
    prior = pseudo * prior;
    
end