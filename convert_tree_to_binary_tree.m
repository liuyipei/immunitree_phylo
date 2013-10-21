function T = convert_tree_to_binary_tree(T,t)

r = find(T == 0);
assert(length(r) == 1);
N = length(t);

T = [t; T]+N;
T(N+r) = 0;
M = size(T,1);
T = process_node(N+r, T);
T = clean_tree(T);
end

function T = process_node(v, T)
    children = find(T == v);
    if length(children) == 1
        % erase self
        T(children) = T(v);
        T(v) = -1;
        T = process_node(children, T);
    elseif length(children) == 2
        T = process_node(children(1), T);
        T = process_node(children(2), T);
    elseif length(children) > 2
        children = children(randperm(length(children)));
        % convert to binary by adding internal nodes
        T = [T; v];
        T(children(2:end)) = size(T,1);
        T = process_node(size(T,1), T);
        T = process_node(children(1), T);
    end
end


function test()
%%
T = [0 1 1 2 3 2 3 2 3 2 3 7 12 12]';
t = [4 6 8 10 5 13 14 9 11]';

%             1
%       2               3
%    4 6 8 10       5  7   9  11
%                     12
%                    13 14
convert_tree_to_newick(convert_tree_to_binary_tree(T,t), true)

%%
T = [0 1 1 2 2]';
% nodes 2,3,4 are leaves
%t = [2*ones(1,10) 3*ones(1,10) 4*ones(1,10)]';
t = [ones(1,10) 3*ones(1,10) 4*ones(1,10) 5*ones(1,10)]';
convert_tree_to_newick(convert_tree_to_binary_tree(T,t), true)
a.T = T;
a.t = t;

T = [0 1 2 2 1]';
% nodes 2,3,4 are leaves
%t = [2*ones(1,10) 3*ones(1,10) 4*ones(1,10)]';
t = [ones(1,10) 3*ones(1,10) 4*ones(1,10) 5*ones(1,10)]';
b.T = T;
b.t = t;

[~, junk]=compare_trees(a.T, a.t, b.T, b.t)

end