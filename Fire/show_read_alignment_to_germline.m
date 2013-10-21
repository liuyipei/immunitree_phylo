function f = show_read_alignment_to_germline(clone, G, src, t)
    ret_gca = -1;
    if ~isempty(src)
        if isnumeric(t)
            [~, o] = sortrows([src t]);
        else
            [~, o] = sortrows([src]);
        end
        src = src(o);
        t = t(o);
    else
        assert(isnumeric(t));
        [~, o] = sortrows(t);        
        t = t(o);
    end
            
    clone = clone(o,:);
    mut = view_read_alignment(clone, G.seq);
    f = figure; 
    X = [G.seq; mut];
    if ~isempty(src)
        X = [X [5 5 5; repmat(5.5-src, 1, 3)]];
    end
    fprintf('image of germline alignment computed.\n');
    h = imagesc(X);  
    fprintf('image of germline alignment shown -- setting labels\n');
    if size(clone, 1) < 10000 % on 45k, matlab encountersseg fault
        set(gca, 'YTick', 2:size(t,1)+1);
        set(gca, 'YTickLabel', t);
        if size(clone,1) > 50
            set(gca, 'FontSize', 8);
        end
        %fprintf('read alignment to germline figure parameters set.\n');
    else
        %fprintf('%d reads : skipping y-tick figure parameters.\n', size(clone, 1));
    end
    
    title('read alignment to germline');
    ret_gca= gca;
end