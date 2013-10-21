function ix = get_reads(data, subject_num, clean, FR, read_frame, vdj, code)

    if ~exist('clean', 'var'), clean = true; end
    if ~exist('subject_num', 'var'), subject_num = []; end
    if ~exist('vdj', 'var'), vdj = zeros(1,3); end
    if ~exist('read_frame', 'var'), read_frame = []; end
    if ~exist('code', 'var'), code = 'ihmmune'; end
    if ~exist('FR', 'var'), FR = []; end
    
    ix = true(length(data.subject_num),1);

    if ~isempty(subject_num)
%        ix = ix & (data.subject_num == subject_num);
        ix = ix & ismember(data.subject_num,subject_num);
    end
    
    if ~isempty(clean) % reads that had no "problems" when aligned to V,J
        ix = ix & (data.trimmed_VJ(:,4) == 0);
    end
    
    if ~isempty(FR)
        ix = ix & (data.FR == FR);
    end
    
    if ~isempty(read_frame)
        iy = mod(sum(data.trimmed_VJ,2),3) == 1;
        if read_frame
            ix = ix & iy;
        else
            ix = ix & ~iy;
        end
    end
    
    if vdj(1) ~= 0
        iv = [1 4:size(data.(code),2)];
        iv = sum(data.(code)(:,iv) == vdj(1), 2) > 0;
        ix = ix & iv;            
    end        
    
    for i=2:3
        if vdj(i) ~= 0
            ix = ix & (data.(code)(:,i) == vdj(i));
        end            
    end    
    
    ix = find(ix);
    
%     % labels are based on subjec_num + FR.
%     [~, ~, labels] = unique([data.subject_num(ix) data.FR(ix)], 'rows');
    
        

end


%     if ~strcmp(code, 'ihmmune')
%         for i=1:3
%             if vdj(i) ~= 0
%                 ix = ix & (data.(code)(:,i) == vdj(i));
%             end
%         end
%     else % ihmmune
