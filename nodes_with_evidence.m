function [evidence evidence_quantity] = nodes_with_evidence(tree,t)
% find subtree with evidence
T = zeros(size(tree,1),2);
T(:,1) = tree(:,1);
T(T(:,1)==-1,1)=1;  % so erased nodes are not a problem
T(t,2) = 1;
T = fill_counts(T);
evidence = T(:,3)>0;  % nodes tied with evidence
if nargin > 1
    evidence_quantity = T(:,3); % the raw number of counts
end
end