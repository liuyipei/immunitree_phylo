function adj_rand = confusion_matrix(real_labels, my_labels)

max_mine = max(my_labels);
my_labels(my_labels == 0) = max_mine + 1;

max_real = max(real_labels);
real_labels(real_labels == 0) = max_real+1;

hitmap = zeros(max_mine+1,max_real+1);
for i=1:length(my_labels)
    hitmap(my_labels(i), real_labels(i)) = hitmap(my_labels(i), real_labels(i)) + 1;    
end
%%
%W = hitmap*diag(1./sum(hitmap)); % normalize
W = hitmap;
[junk I] = max(W); % for every true sequence s, I(s) is where most its read reside.
[B I_] = unique(I, 'first');  % I_ = the indexes of a unique (and sorted) set of elements in I
sigma = I(sort(I_)); % the unique set of elements in I, maintaining in the order they appear in I
not_in_sigma = setdiff( 1:(max_mine+1), sigma);
sigma = [sigma not_in_sigma];
if nargout == 0
    image(W(sigma, :)); 
    set(gca, 'YTick', [1:length(sigma)]);
    set(gca, 'YTickLabel', num2cell(sigma));
    xlabel('real labels');
    ylabel('my labels');
end
%%  calculating the adjusted RAND index
W = [W sum(W,2); sum(W,1) sum(sum(W))];
W_1 = W-1;
Z = 0.5*W.*W_1;

% a, the number of pairs of elements in S that are in the same set in X and in the same set in Y
% b, the number of pairs of elements in S that are in different sets in X and in different sets in Y
% c, the number of pairs of elements in S that are in the same set in X and in different sets in Y
% d, the number of pairs of elements in S that are in different sets in X and in the same set in Y

a = sum(sum(Z(1:end-1, 1:end-1)));
c = sum(Z(end, 1:end-1)) - a;
d = sum(Z(1:end-1, end)) - a;
b = Z(end,end)-a-d-c;
adj_rand = 2*(a*b-c*d)/( (a+d)*(d+b) + (a+c)*(c+b) );
if nargout == 0,
    title(sprintf('rand index is %.3f', adj_rand));
end
end