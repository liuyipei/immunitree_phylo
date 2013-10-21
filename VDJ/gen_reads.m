function [reads t nodes sequences T] = gen_reads(variants, N, noise, indel, flat)
% Generates N reads total from variants.
% returns reads and assignment of reds to variants (t).
% assignment to variants is sampled using DIR(2*alpha/K)
% where K is the number of variants.

if nargin <5, flat = 1; end

M = length(variants);
alpha = 2*M;

if isempty(variants) || isempty(variants{1})
    reads = []; return;
end

if variants{1}(1) > 5
    map('ACGTNacgtn') = [1:5 1:5];
    for i=1:M
        variants{i} = map(variants{i});
    end
end
    
[t counter] = partition_reads(N, M, alpha);
assert(length(t) == N);
t = sort(t);

reads = cell(N,1);
labels = cell(N,1);

if flat 
    for i=1:N
        R = gen_one(variants{t(i)}, noise, indel);
        reads{i} = R;
        labels{i} = sprintf('id %d source %d', i, t(i));
    end
else    
    % generate from DPtree
    me = 0.5; 
    alpha = 0.5;
    sigma = 0.025;
    epsilon = 0;
    ix = [0 cumsum(counter)];
    for s=1:length(counter)
        [~, sequences{s}, T{s}, nodes(ix(s)+1:ix(s+1))] = ...
            generate_data(counter(s), variants{s}, me, alpha, epsilon, sigma);        
        fprintf('-------------------\nClone number %d:\n\n', s);
        view_tree(T{s});
    end
    for i=1:N
        R = gen_one(sequences{t(i)}(nodes(i),:), noise, indel);
        reads{i} = R;
        labels{i} = sprintf('id %d source %d', i, t(i));
    end
end

end


function [c, counter] = partition_reads(N, M, alpha)
    counter = zeros(1,M);
    c = zeros(N,1);
    for i=1:N
        r = mnrnd(counter+(1/M)*alpha);
        c(i) = r;
        counter(r) = counter(r)+1;
    end
end


function [V change] = gen_one(X, noise, indel)
    V = X;
    assert(all(V>0));

    % deletions
    ix = find(rand(1,length(V))<1-indel);    
    V = V(ix);
    change = ~isempty(ix);
    
    %insertions
    ix = find(rand(1,length(V))<indel);        
    V = [V; zeros(1,length(V))];
    V(2*ix) = ceil(4*rand(1,length(ix)));       
    V = V(V>0)';
    change = change || ~isempty(ix);
    assert(all(V>0));
    
    % noise
    ix = find(rand(1,length(V))<noise*4/3);
    V(ix) = ceil(4*rand(1,length(ix)));       
    change = change || ~isempty(ix);
    assert(all(V>0));
end

function r = mnrnd(p)
    p = p/sum(p);
    P = cumsum(p);
    z = rand;
    r = find(z < P, 1, 'first');
end