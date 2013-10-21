function [clone germline id src rep st h] =  play_with_clone(clone_file, st, data)
    
    % load clone .fa file    
    if ischar(clone_file)
        W = fastaread(clone_file);    
    else 
        W = clone_file;
    end
    
    W_ = {W(2:end).Header};    
    W_ = cell2mat(cellfun(@(x) [x '_'], W_, 'uniformoutput', false));
    [str_id, id, src rep t] = strread(W_,'%s%d%*s%d%*s%d%*s%d','delimiter','_');
           
    clone = cell2mat({W(2:end).Sequence}');
    map('ACGTN') = 1:5;
    clone = map(clone);

    % create germline
    num_seperators = length(find(W(1).Header == '_'));
    assert(num_seperators == 13);    
    germline = W(1);
    germline.seq = map(germline.Sequence);
    eaten = zeros(1,6);
    [d eaten(1) eaten(2) eaten(3) eaten(4) eaten(5) eaten(6) lastV firstJ n1 n2] = ...
        strread([W(1).Header '_'],'%*s%*s%*s%s%d%d%d%d%d%d%d%d%s%s','delimiter','_');    
    n1 = n1{1}; n2 = n2{1}; d=d{1};

    germline.full = germline.seq;
    germline.full(lastV+(1:length(n1))) = map(n1);
    germline.full(firstJ-(length(n2):-1:1)) = map(n2);
    
    germline.eaten = eaten;
    germline.d = d;
    germline.n1 = n1;
    germline.n2 = n2;
    germline.lastV_firstJ = [lastV firstJ];
    germline.frame_shift = mod(eaten(1)+eaten(6) + length(germline.seq) - 1, 3);

    % load clone .mat file
    if nargin<2 && (nargout > 5 || nargout == 0)
        assert(isstr(clone_file));
        mat_file_name = [clone_file(1:end-3) '.mat'];    
        Z = load(mat_file_name);
        st = Z.a;        
    end

    % if st is not empty, organize the nodes of the tree according to DFS
    % order
    if exist('st', 'var')
        st = order_nodes(st);
    end    
    
    
    % otherwise assumes that data is the t5.fa file *after* PCR
    if exist('data', 'var') && ~isempty(data)
        id = cell2mat(str_id);
        [~,id] = ismember(id, data.read_id, 'rows');
    end
%    keyboard;

if nargout == 0 || nargout == 7
    
    assert(exist('st', 'var') ~= 0);
    if ischar(clone_file)
        filename = clone_file(find(clone_file == '/', 1, 'last')+1:end);
    else 
        filename = '';
    end
    
    str = sprintf('  %s\nN1/2: [%d %d]; ', filename, length(n1), length(n2));
    eaten_str = sprintf(' %d', eaten(2:5)); 
    str = [str sprintf(' eaten: [%s]', eaten_str(2:end))];
    str(str == '_') = ' ';

    h = visualize_tree(st.tree, st.sequences, [st.t src], clone, 4, germline);  
    print_str = false;
    if ~isempty(h) && print_str
        h.Nodes(1).Label = [h.Nodes(1).Label str];
    end
    show_silent(clone, germline, src);
    show_read_alignment_to_germline(clone, germline, src, st.t);    
    set(gcf, 'Position', [6 281 1126 246]);
end

end





function show_silent(clone, G, src)
    SM = silent_map();
    
    mut = view_read_alignment(clone, G.seq);
    mut(G.seq == 5) = 5; % ignore mutations from junction regions.
    
    if G.frame_shift ~= 0
        fprintf('can''t show silent mutation because clone is out of frame\n');
    end
    
    codons = seqs2codons(clone(1,1:end-mod(size(clone,2),3)));
    codons = reshape(repmat(codons, 3,1), 1, []);    
    
    % for every (location, letter,source)
    H = zeros(size(mut,2), 5, 4); % sites x letters x sources
    G = zeros(size(mut,2), 5);
    for j=1:length(codons)
        G(j,1:4) = double(SM(codons(j), mod(j-1,3)+1, :));        
    end
    
    for j=1:size(mut,1)
        ix = find(mut(j,:) ~= 5);
        H(sub2ind(size(H), ix, mut(j,ix), src(j)*ones(1,length(ix)))) = ...
            H(sub2ind(size(H), ix, mut(j,ix), src(j)*ones(1,length(ix)))) + ...
            G(sub2ind(size(G), ix, mut(j,ix)));
    end
    
    ix = max(sum(H  ,3),[],2) < size(mut,1);
    H_ = max(sum(H>0,3),[],2);
    iy = H_>=2;
    fprintf('Silent Mutations:\n');
    [find(ix&iy)'; H_(ix&iy)']
end



