function clone_dist = calculate_average_read_distance(T, nMuts, t, src)
    h = T(:,2);
    assert(h(1) == 0);
    if isempty(t), h(h>0) = 1; end  % count each non-empty node as 1 read.
    nMuts(2) = 1;
    D = calculate_tree_distances(T, nMuts);
    B = calc_frequency_of_clone_pairs(h);

% alternative to nMuts(2) = 1 
%     D = D(2:end, 2:end);
%     B = B(2:end, 2:end);
    clone_dist = B(:)'*D(:);

    if nargin == 4
        clone_dist = [zeros(1,4) clone_dist];
        for s=1:4
            ix = (src == s);
            h = hist(t(ix), 1:size(T,1))';
            B = calc_frequency_of_clone_pairs(h);
            clone_dist(s) = B(:)'*D(:);
        end
    end
end


function D = calculate_tree_distances(T, nMuts)
    M = size(T,1);
    DG = sparse(1+T(:,1), 2:M+1, nMuts, M+1, M+1);
    DG = DG(2:end,2:end);
    DG = DG+DG';
    D = graphallshortestpaths(DG, 'directed', false);

end


function B = calc_frequency_of_clone_pairs(h)
    M = length(h);
    %if sum(h) == 1, B = 0; return;  % do we want to return NaN or 0 if
    %there is only one read?
    [I,J] = meshgrid(1:M, 1:M);
    B = h*h';
    B(I<J) = 0;
    B(I == J) = h.*(h-1) / 2;
    assert(sum(B(:)) == sum(h)*(sum(h)-1)/2);
    B = B/sum(B(:));
end