function P = Ztest(H, W)

% not sure why, but this is how the website does it at 
% http://people.chgv.Isrc.duke.edu/~dg48//metap.php
%W = W-2;  % you can remove if necessary.

% get the Z scores
Z = norminv(H);
W(isnan(Z)) = 0;
Z(isnan(Z)) = 0;
Z(isinf(Z)) = 3; % arbitrary number?

Z = Z .* W;
K = sqrt(sum(W.^2,3));
X = sum(Z,3)./K;
% X ~ N(0,1)
P = normcdf(X);


end


%%