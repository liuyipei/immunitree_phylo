function [Q pi total_mutations] = estimate_mutation_parameters(tree, sequences, Q, rate_class)

assert(sum(tree(:,1)==-1) == 0);
assert(tree(1,1) == 0);
nClasses = size(Q,3); % this is the only reason we really need Q
L = size(sequences,2);
B = size(Q,1);

epsilon = 0.001;
pseudo = 1000;
pseudo = pseudo*(epsilon+(1-B*epsilon)*eye(B));
pseudo = pseudo(:,:,ones(1,nClasses));
counts = zeros(size(Q));


for l=1:L            
    ix = (sequences(2:end, l) ~= sequences(tree(2:end,1),l));
    for u=(1+find(ix))'
        par = tree(u,1);
        a = sequences(par,l);
        b = sequences(u,l);
        % remember, Q is indexed (child,parent,class)
        counts(b,a,rate_class(l)) = counts(b,a,rate_class(l))+1;
    end
    for b=1:B
        counts(b,b,rate_class(l)) = counts(b,b,rate_class(l)) + sum(~ix & sequences(2:end,l) == b);
    end
end

Q = pseudo + counts;
Q = reshape(Q,B,[]);

sample = true;
if sample
    % sample Q
    for u=1:size(Q,2)
        Q(:,u) = drchrnd(Q(:,u)', 1);
    end
else
    % find the MAP Q
    sumQ = sum(Q,1);
    Q = Q./sumQ(ones(1,4),:);  % normalize columns
end

Q = reshape(Q,B,B,nClasses);

% mutation stat
counts = sum(counts, 3);
total_mutations = sum(counts(:))-sum(diag(counts));

% compute pi (prior over root)
pseudo_pi =10;
pi = hist(sequences(1,:), 1:B) + pseudo_pi;
pi = pi/sum(pi);

end



% 