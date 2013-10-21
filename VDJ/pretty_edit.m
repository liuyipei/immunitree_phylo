function str = pretty_edit(edit)    
    str = '';
    if length(edit) == 0
        return;
    end
    if edit(1) > 5  % this is already the formatted string 
        str = edit; 
        return; 
    end
    
    map = ['D', 'M', 'I', 'E'];
    curr_ch = edit(1);
    count = 0;
    for i=1:length(edit)
        if edit(i) ~= curr_ch
            str = [str int2str(count) map(curr_ch+2) '-'];
            count = 0;
            curr_ch = edit(i);
        end
        count = count+1;
    end
    str = [str int2str(count) map(curr_ch+2)];
end
