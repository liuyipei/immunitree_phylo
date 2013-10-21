function prior = generate_sticky_prior(epsilon, B, pseudo)
    prior = pseudo * (epsilon + (1-B*epsilon)*eye(B));
end