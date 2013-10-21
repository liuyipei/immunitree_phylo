function visualize(c_, aligned, edt, T_, seqs_, t_, time_point)
    
    g = figure('WindowButtonDownFcn', @debug_me, 'DeleteFcn', @close_me);
    N = length(c_);
    M = max(c_);
    dict = 'ACGTN';
    
    sum_edt = zeros(N,1);
    for i=1:N, sum_edt(i) = sum(edt{i} ~= 0); end
    
    for k=1:min(16, M)
        h = subplot(4,4,k);
        ix = find(c_==k);
        aligned_k = cell2mat(aligned(ix));
        view_read_alignment(aligned_k);
        title(sprintf('clone %d (length %db)', k, size(aligned_k, 2)));
        xlim([0 size(aligned_k, 2)]);
        ylim([0 length(ix)]);
        ax = axis;
        text(0.85*ax(2),0.6*ax(4),sprintf('avg\nedit\ndist\n%.1f', mean(sum_edt(ix))));
        set(h, 'UserData', k);
    end
    subplot(4,4,2);
    tmp = ylim;
    text(-10,1.25*tmp(2),'Number of reads possesing each allele, on each location')
    legend({'rarest', '  : ', '  : ', 'common', 'missing'})


function close_me(src, evt)
    h = get(src, 'Userdata');
    if isempty(h)
        h = zeros(2,1);
    end
    if h(1) ~= 0
        close(h(1));
    end    
    if h(2) ~= 0
        close(h(2));
    end
end
    
    
function debug_me(src, evt)
    h = get(src, 'Userdata');
    if isempty(h)
        h = zeros(2,1);
    end

% show sequence    
    pos1 = [1500 1000 1000 400];
    if h(1) ~= 0
        pos1 = get(h(1), 'position');
        close(h(1));
    end
    
    o  = get(gca, 'CurrentPoint');
    i  = get(gca, 'Userdata');
    ix = find(c_==i);
    al = cell2mat(aligned(ix));
    st = max(1,round(o(1,1))-20);
    en = min(size(al, 2), round(o(1,1))+20);
    [~,h(1)] = seqlogo(dict(al), 'startat', st,'endat', en);
    set(h(1), 'Position', pos1);
    set(h(1), 'DeleteFcn', {@del_me, 1});
    
% Show tree    
    if strcmp(get(src, 'SelectionType'), 'extend')
        pos2 = [];
        if h(2) ~= 0
            pos2 = get(h(2), 'position');
            close(h(2));
        end
        time_stamp_data = [t_{i}; time_point(ix)]';
        j = view_tree(T_{i}, seqs_{i}, 1, time_stamp_data);                
        h(2) = gcf;
        if ~isempty(pos2)
            set(h(2), 'Position', pos2);
        end    
        set(h(2), 'DeleteFcn', {@del_me, 2});
                
    end
    
    set(src, 'Userdata', h);
    
    function del_me(a,b,c)
        junk = [1 1]'; junk(c) = 0;
        set(src, 'UserData', junk .* get(src, 'UserData'));            
    end
end

end