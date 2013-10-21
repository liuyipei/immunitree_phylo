function print_text_on_image(P)
    M = size(P,1);
    assert(M == size(P,2));
    
    for l=1:M
        for l_ = l+1:M
            vals = [P(l,l_) P(l_,l)];
            vals = sort(vals);
            str = '';
            if vals(2) ~= 0
                str = sprintf('%d',vals(2));
            end            
            if vals(1) ~= 0
                str = [str sprintf('|%d',vals(1))];
            end
            if ~isempty(str) 
                h = text(l_-0.2,l, str); 
%                set(h, 'color', [1 0 1]);
            end
            
        end
    end    
end