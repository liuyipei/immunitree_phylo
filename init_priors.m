%function priors = init_priors(priors_)
% Sets the priors according to the values in priors_ and fills the missing
% values with defaults.  
% Explanation of structure fields:
%   is_codon: true/false if working in codon/nucleotide mode
%   R: Dirichlet weights on the read noise 4x4 transition matrix
%   NT: Dirichlet weights on the mutation nucleotide 4x4 transition factor
%   AA: Dirichlet weights on the mutation amino-acid 4x4 transition factor
%       (only applicable in codon mode)
%   pi: Dirichlet weights on the prior over the root sequence
%   decay: gaussian prior over the decay parameter in mutation model
%   br_shape: birth rate shape parameter (scale parameter is 1)
%   dr_shape: death rate shape parameter (scale parameter is 1)
function priors = init_priors(priors_)

    % first initialize the prior stuct to default values
    priors.is_codon = true;
    priors.R = generate_sticky_prior(0.0001, 4, 630000);
    priors.NT = generate_sticky_prior(0.001, 4, 630);
    priors.AA = generate_sticky_prior(1/21, 21, 210);
    priors.pi = 100*ones(1,64);    
    priors.decay = [0 1e-3];  % mean and std
    priors.nClasses = 3;
    priors.br_shape = 1;
    priors.dr_shape = 0.75;

    if isfield(priors_, 'is_codon') && priors_.is_codon == false
        % default values for nucleotide mode
        priors.pi = 100*ones(1,4);    
        priors.is_codon = false;
    end
    
    if isfield(priors_, 'is_decay') && priors_.is_decay == false
        % default values for no decay
        priors.decay = [0 0];  % mean and std
    end

    % override default values with user input values (in 'priors_')
    priors = struct_override(priors, priors_);
end