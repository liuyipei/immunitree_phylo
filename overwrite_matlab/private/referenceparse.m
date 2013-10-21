function [data,gptext,ln] = referenceparse(data,gptext,ln,record_count)
%REFERENCEPARSE parses the reference entries of GenBank files

% Copyright 2003-2006 The MathWorks, Inc.
% $Revision: 1.1.12.3 $   $Date: 2006/06/16 20:07:01 $

ref_count=1;
data(record_count).Reference = [];
while ~matchstart(gptext(ln,:),'FEATURES') && ~matchstart(gptext(ln,:),'COMMENT') && ~matchstart(gptext(ln,:),'ORIGIN')
    
    data(record_count).Reference{ref_count}.Number = '';
    data(record_count).Reference{ref_count}.Authors = '';
    data(record_count).Reference{ref_count}.Consrtm = '';
    data(record_count).Reference{ref_count}.Title = '';
    data(record_count).Reference{ref_count}.Journal = '';
    data(record_count).Reference{ref_count}.MedLine = '';
    data(record_count).Reference{ref_count}.PubMed = '';
    data(record_count).Reference{ref_count}.Remark = '';
    

    %REFERENCE - Mandatory
    [s,f,t] = regexp(gptext(ln,:),'REFERENCE\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Number = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end

    %AUTHORS - Mandatory -- though we found examples where this was missing
    [s,f,t] = regexp(gptext(ln,:),'AUTHORS\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Authors = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end
    
    while ~matchstart(gptext(ln,:),'TITLE') && ~matchstart(gptext(ln,:),'JOURNAL') && ~matchstart(gptext(ln,:),'CONSRTM')
        data(record_count).Reference{ref_count}.Authors=strvcat(data(record_count).Reference{ref_count}.Authors, strtrim(gptext(ln,:))); %#ok
        ln=ln+1;
    end
    
    %CONSRTM - Optional    
    [s,f,t] = regexp(gptext(ln,:),'CONSRTM\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Consrtm = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end
    while ~matchstart(gptext(ln,:),'TITLE') && ~matchstart(gptext(ln,:),'JOURNAL')
        data(record_count).Reference{ref_count}.Consrtm=strvcat(data(record_count).Reference{ref_count}.Consrtm, strtrim(gptext(ln,:))); %#ok
        ln=ln+1;
    end

    %TITLE - Optional    
    [s,f,t] = regexp(gptext(ln,:),'TITLE\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Title = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end
    while ~matchstart(gptext(ln,:),'JOURNAL')
        data(record_count).Reference{ref_count}.Title=strvcat(data(record_count).Reference{ref_count}.Title, strtrim(gptext(ln,:))); %#ok
        ln=ln+1;
    end

    %JOURNAL - Mandatory
    [s,f,t] = regexp(gptext(ln,:),'JOURNAL\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Journal = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end
    
    if matchstart(gptext(ln,:),'REFERENCE')
        ref_count=ref_count+1;
        % next reference
        continue
    end

    if matchstart(gptext(ln,:),'COMMENT') || matchstart(gptext(ln,:),'FEATURES') || matchstart(gptext(ln,:),'BASE COUNT')
        % done with references
        break
    end

    while ~matchstart(gptext(ln,:),'MEDLINE') && ~matchstart(gptext(ln,:),'PUBMED') && ~matchstart(gptext(ln,:),'REMARK')
        data(record_count).Reference{ref_count}.Journal=strvcat(data(record_count).Reference{ref_count}.Journal, strtrim(gptext(ln,:))); %#ok
        ln=ln+1;
        if matchstart(gptext(ln,:),'REFERENCE') || matchstart(gptext(ln,:),'COMMENT') || matchstart(gptext(ln,:),'FEATURES') || matchstart(gptext(ln,:),'BASE COUNT')
            % done with references
            break
        end
    end

    %MEDLINE - Optional    
    [s,f,t] = regexp(gptext(ln,:),'MEDLINE\s+(\d)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.MedLine = gptext(ln,t{1}(1):t{1}(2));
        ln=ln+1;        
    end

    %PUBMED - Optional    
    [s,f,t] = regexp(gptext(ln,:),'PUBMED\s+(\d+)');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.PubMed = gptext(ln,t{1}(1):t{1}(2));
        ln=ln+1;
    end

    %REMARK - Optional    
    [s,f,t] = regexp(gptext(ln,:),'REMARK\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Remark = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end
    while ~isempty(s) && ~matchstart(gptext(ln,:),'COMMENT') && ~matchstart(gptext(ln,:),'FEATURES') && ~matchstart(gptext(ln,:),'BASE COUNT') && ~matchstart(gptext(ln,:),'REFERENCE')
        data(record_count).Reference{ref_count}.Journal=strvcat(data(record_count).Reference{ref_count}.Journal, strtrim(gptext(ln,t{1}(1):t{1}(2)))); %#ok
        ln=ln+1;
        if matchstart(gptext(ln,:),'COMMENT') || matchstart(gptext(ln,:),'FEATURES') || matchstart(gptext(ln,:),'BASE COUNT'), break, end
    end

    % next reference
    ref_count=ref_count+1;
end