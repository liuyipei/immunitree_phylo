function Q = generate_stochastic_matrix(prior, nClasses, use_gamma)    
    if ~exist('use_gamma', 'var'), use_gamma = true; end
    
    Q = randg(prior(:,:,ones(1,nClasses)));     % result is unnormalized!
        
    if ~use_gamma   % Dirichlet.  Normalize columns
        sumQ = sum(Q,1);
        Q = Q./sumQ(ones(1,size(Q,1)), :);
    end
end