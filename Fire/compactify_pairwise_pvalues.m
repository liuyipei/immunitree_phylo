function Q = compactify_pairwise_pvalues(P)
% Show a restricted view on P
Q = cell(length(P),1);
for j=1:length(P)
    Q{j} = compactify_pairwise_pvalues_inner(P{j});
end
 

    
end

function Q = compactify_pairwise_pvalues_inner(P)
    [codes, ~, doubles] = get_tuple_codes();
    Q = ones(4);
    for i=1:length(doubles)
        src = find(codes(doubles(i),:));
        Q(src(1), src(2)) = min(P(src(1), 4+i), P(4+i, src(1))); 
        Q(src(2), src(1)) = min(P(src(2), 4+i), P(4+i, src(2)));
    end
end