function block = getfileblock(filename,range,delimiter)
% GETFILEBLOCK reads a block of a file
%
%   GETFILEBLOCK(FILENAME,RANGE,DELIMITER) reads a block of text from file
%   FILENAME. The function assumes that the file is made up of entries with
%   a delimiter DELIMITER marking the start or end of an entry. The block
%   will start at RANGE(1) and stop at RANGE(2).
%

%   Copyright 2002-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2010/12/22 16:19:19 $%

% set start and stop and make sure they are reasonable
start = floor(range(1));
if start < 1
    start = 1;
end
if numel(range) > 1
    stop = ceil(range(2));
else
    stop = start;
end

fid = fopen(filename,'r');
% read in blocks of 16MB at a time
blockSize = 2^24;

hNdx = [];
pos = 0;
count = blockSize;
bufferLen = length(delimiter)-1;
buffer = '';
while count == blockSize
    [str count] = fread(fid,blockSize,'*char');
    str = str';
    if length(str) < bufferLen
        break
    end
    found = strfind([buffer str(1:bufferLen)],delimiter);
    if ~isempty(found)
        hNdx = [hNdx (pos-bufferLen + found)];
    end
    found = strfind(str,delimiter);
    if ~isempty(found)
        hNdx = [hNdx (pos + found)];
    end
    pos = ftell(fid);
    if numel(hNdx)> stop
        break;
    end
    buffer = str(end-bufferLen+1:end);
end
frewind(fid);


numItems = numel(hNdx);

if start > numItems
    error(message('bioinfo:getfileblock:StartTooBig', numItems));
end
stop = min(numItems,stop);

fclose(fid);
if numItems == 0
    block = '';
else

    fid = fopen(filename,'r');
    chunkStart = hNdx(start)-1;
    if isfinite(stop) && stop < numItems
        chunkEnd = hNdx(stop+1)-1;
        chunkSize = chunkEnd - chunkStart;
    else
        chunkSize = Inf;
    end
    fseek(fid,chunkStart,-1);
    try
        block = fread(fid,chunkSize,'*char')';
    catch theErr
        if strcmpi(theErr.identifier,'MATLAB:nomem')
            error(message('bioinfo:getfileblock:BlockTooBig'));
        else
            rethrow(theErr);
        end
    end
    fclose(fid);
end

