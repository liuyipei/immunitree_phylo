function [v j big_hist] = collect_VJs(data, ix, method, thresh)

%    ix = get_reads(data, cur_sample);

    % V-J distribution
    spectrum = max(data.(method));
    spectrum(2) = spectrum(2)+1;
    big_hist = zeros(spectrum(1:3));
    iv = [1 4:size(data.(method),2)];
               
%     if strcmp(method, 'ihmmune')
%         ihmmune = true;
%         erased = zeros(1,max(spectrum));
%         % These Vs were marked redundant by "find_redundant_Vs.m" 
%         % because there is another V that always appears with them
%         redundant = [6 7 24 26 34 58 65 70 71 74 109 116 135 138 140 146 ...
%             150 152 153 170 174 176 178 190 204 219 221 224 227 253 305 ...
%             311 315 322 324];
%         erased(redundant) = 1;
%     else
%          ihmmune = false;
%     end
%         
    
    vdj = data.(method);
    for i=ix'        
        v = vdj(i,1); d = vdj(i,2); j = vdj(i,3);
        if v == 0 || j == 0, continue; end
        if d == 0, d = spectrum(2); end
        vs = vdj(i, iv);
        for v=vs
            if v == 0, break; end
%             if ihmmune && (erased(v) == 1), continue; end
            big_hist(v,d,j) = big_hist(v,d,j) +1;
        end
    end
    
    big_hist = squeeze(sum(big_hist,2));
    if nargout == 0, imagesc(big_hist'); end
    [sorted_hist, ordered_VJs] = sort(big_hist(:), 'descend');
    if thresh >=0
        nTop = find(sorted_hist < thresh,1)-1;  % the last included entry
    else
        nTop = min(-thresh, length(sorted_hist));
    end
    if isempty(nTop), nTop = length(sorted_hist); end

    [v j] = ind2sub(size(big_hist), ordered_VJs(1:nTop));
%    fprintf('Found %d VJ-clones with at least %d reads\n', nTop-1, thresh);
end


