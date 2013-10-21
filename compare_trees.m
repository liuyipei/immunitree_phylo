function [soft jaccard]= compare_trees(T1, t1, T2, t2)

N = length(t1);
assert(N == length(t2));

% for each edge in T1, construct a hamming vector (or two) representing the
% splits.  do the same for each edge in T2. 
B1 = get_splits(T1(:,1),t1);
M = size(B1,1);
B1 = [B1; ~B1];
B2 = get_splits(T2(:,1),t2);

% How many vectors of T2 are in T1?
B = intersect(B1,B2, 'rows');
D = 1-(double(B1)*double(B2)' + double(~B1)*double(~B2)')/N;
D = [D(1:M,:) D(M+1:end,:)]; 
jaccard = size(B,1)/(M+size(B2,1)-size(B,1));
soft = mean(min(D,[],2));
end

function order = get_topological_order(T)
    M = size(T,1);
    if M == 1,  order = 1; return; end
    DG = sparse(1+T, 2:M+1, true, M+1, M+1);
    order = graphtopoorder(DG(2:end, 2:end));
end

function B = get_splits(T,t);
    order = get_topological_order(T);
    rev_order = fliplr(order);
    B = false(length(T),length(t));
    for v=rev_order(1:end-1)
        B(v,(t==v)) = true;
        B(T(v),B(v,:)) = true;
    end
    % add the reads at the root 
    B(rev_order(end),:) = true;
    B = unique(B,'rows');
end

% map = 'ACGT';
% Z = []; for i=1:size(reads,1), Z(i).Sequence = map(reads(i,:)); Z(i).Header = sprintf('id_%d_node_id_%d\n', i, b.t(i)); end
% fastawrite('tandy.fa', Z);
% fasta2relaxedPhylip.pl tandy.fa
% raxmlHPC -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s ex_al -n TEST
% FastTree -gtr -nt alignment_file > tree_file 
