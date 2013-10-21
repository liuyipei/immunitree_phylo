function [rates ll] = MH_rates(tree, F, rates, priors, action)

global greedy

if nargin <4, action = 1; end

% DONE: work only on the part of the tree with no evidence attached
% I actually can't do that!  Because I need to condition on the fact that
% the part of the tree with no evidence has the given number of "live"
% cells (or compute the distribution over the number of live cells).
% Commenting out for now.

% evidence = nodes_with_evidence(tree,t);
% tree(~evidence,1) = -1;
% tree = clean_tree(tree);

% prior on rates (gamma with scale=1 and the following shapes):
br_shape = priors.br_shape; %=1;
dr_shape = priors.dr_shape; %=0.75;

nBirths = size(tree,1)-1;
nDeaths = sum(tree(:,3) < F);
total_life_time = sum(tree(:,3)-tree(:,2));

% change birth-rate
if action ~= 0
    if greedy == -1 || action == -1
        rates.birth = (nBirths+br_shape)/(1+total_life_time);
        rates.death = (nDeaths+dr_shape)/(1+total_life_time);
    else
        rates.birth = randg(nBirths+br_shape)/(1+total_life_time);
        rates.death = randg(nDeaths+dr_shape)/(1+total_life_time);
    end
end

ll = nBirths*log(rates.birth) + nDeaths*log(rates.death) - (rates.death+rates.birth)*total_life_time;

end

