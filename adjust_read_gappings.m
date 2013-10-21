% input:
% T is Mx1 array.  T(i) is the parent of node i.
% singles is MxBxL, where L is the number of sites (length of each seq).
%    At the end this will store the marginal probabilities for each node of
%    the tree.  Now it stores some prior information on each node.
% pairwise: a BxBxL array with pairwise 
%    potentials for each edge in the tree. pairwise(i,j,l) is the probability
%    a parent with base j generated a child with base i, on location l.
%    exception:  If pairwise is a BxB matrix, then convert to BxBxL by:
%        Q = Q';  %reverse parent_child to child_parent potential
%        Q = Q(:,:,ones(1,L));
function reads_ = ...
    adjust_read_gappings(tree_, sequences_, t_, reads_, raw_from_uniq, uniq_from_raw)

T = tree_(:,1);   
M = length(T);

% Phase one - compute unnormalized conditionals:
% go over tree in reverse topological ordering
DG = sparse(1+T, 2:M+1, true, M+1, M+1);
if M == 1, 
    order = 1;
else
    order = graphtopoorder(DG(2:end, 2:end));
end

adjust_status = struct(...    
    'readgaps_realigned', 0);
for r = raw_from_uniq(:)' % for a representative from each unique sequence
    node_seq = sequences_(t_(r),:); % sequence of the node to which it belongs
    diff = find(reads_(r,:) ~= node_seq);
    for d = 1:(length(diff)-1)
        for e = (d+1):length(diff)
            positions = diff(d):diff(e);
            r_non5_pos = non5(reads_(r,positions));
            s_non5_pos = non5(node_seq(positions));
            if isequal(r_non5_pos, s_non5_pos)
                reads_(r,positions) = node_seq(positions);
                adjust_status.readgaps_realigned = adjust_status.readgaps_realigned + 1;
                break;
            else
                r_s_min_size = min(length(r_non5_pos), length(s_non5_pos));
                if ~isequal(r_non5_pos(1:r_s_min_size), s_non5_pos(1:r_s_min_size))
                    % prefixes do not match, so we can move to the next
                    % START position
                    break;
                end
                % otherwise, we move on to the next END position
            end            
        end
    end
end

for rr = 1:length(t_) % for each generic
    r = raw_from_uniq(uniq_from_raw(rr)); % get the unique representative
    if r == rr
        continue
    else
        reads_(rr,:) = reads_(r,:); % make sure changes are carried over
    end
    
end

end

function non5= non5(x)
    non5= x(x~=5);
end