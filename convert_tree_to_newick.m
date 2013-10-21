function str = convert_tree_to_newick(T, nexus)
    if ~exist('nexus', 'var'), nexus = false; end
    M = size(T,1);
    r = find(T(:,1) == 0);
    assert(length(r) == 1);
    DG = sparse(1+T(:,1), 2:M+1, true, M+1, M+1); 
    DG = DG(2:end, 2:end);        
    str = write_node(r, DG);
    
    if nexus
        fprintf('#NEXUS\n');
        fprintf('BEGIN TREES;\n');
        fprintf('\tTRANSLATE\n');
        for i=1:r-1
            fprintf('\t\t%d\t%d,\n', i, i);
        end
        fprintf('\tUTREE PAUP_1=\n');
        fprintf('%s;\n', str);
        fprintf('ENDBLOCK;\n');          
    end
end

function str = write_node(num, DG)
    children = find(DG(num,:));
    assert(length(children) == 0 || length(children) == 2);
    str = '';
    for t=children
        str=[str write_node(t,DG) ','];
    end
    if ~isempty(str) 
        str = sprintf('(%s)', str(1:end-1)); 
    else % print node_id only for leaves
        str = sprintf('%s%d', str, num);
    end
end


% TODO:  If there is only one child, delete it?
% TODO:  If there are more than two children, then what?