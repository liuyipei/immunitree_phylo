function prior = generate_sticky_with_indel_prior(epsilon, B, pseudo, hm_indel_perc)
    % the last columns and rows are now correct
    hm_indel_rate = (B-1) * epsilon * hm_indel_perc; %(B-1)*epsilon chance of change, then scale by indel chance
    % allocate the indel prob mass away on the right and bottom fringe
    prior = generate_sticky_prior(hm_indel_rate , B, pseudo); 
    % allocate the rest
    prior(1:(B-1), 1:(B-1)) = generate_sticky_prior(...
        epsilon, B-1, pseudo*(1-hm_indel_rate ));
    
end