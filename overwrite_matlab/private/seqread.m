function [data, filetype] = seqread(filename, varargin)
%SEQREAD reads a sequence from a file.
%
%   SEQ = SEQREAD(FILENAME) reads a sequence from file FILENAME, returning
%   the data in the file as a structure. FILENAME can also be a URL or
%   MATLAB character array that contains the text of a supported file format.
%   Supported formats are EMBL, FASTA, GENBANK, GENPEPT, PDB and raw text.
%
%   [SEQ, FORMAT] = SEQREAD(FILENAME) returns the format of the input. 
%
%   See also EMBLREAD, FASTAREAD, GENBANKREAD, GENPEPTREAD, PDBREAD.

%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision: 1.1.8.6 $  $Date: 2010/12/22 16:19:29 $

%  SEQREAD(...,'preambletext', true) gets preamble info for genpept, 
%  genbank, and embl 

% get the text into memory
% in a future version we may accept also cells

getPreambleText = false;

% process input arguments
if  nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:seqread:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'preambletext',''};
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:seqread:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:seqread:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % 'preambletext'
                    getPreambleText = opttf(pval);
            end
        end
    end
end


if ~ischar(filename)
    error(message('bioinfo:seqread:InvalidInput'))
end
dispname = filename;
if size(filename,1)>1  % is padded string

    seqtext = filename;
    seqtext(:,end+1) = sprintf('\n');
    seqtext = seqtext';
    seqtext = seqtext(:)';
    dispname = 'MATLAB String';
    % try then if it is an url
elseif (strfind(filename(1:min(10,end)), '://'))
    if (~usejava('jvm'))
        error(message('bioinfo:seqread:NoJava'))
    end
    try
        seqtext = urlread(filename);
    catch allExceptions %#ok<NASGU>
        error(message('bioinfo:seqread:CannotReadURL', filename));
    end
    % ftext = strread(ftext,'%s','delimiter','\n');

    % try then if it is a valid filename
elseif  (exist(filename,'file') || exist(fullfile(pwd,filename),'file'))
    try
        fid = fopen(filename,'rt');
        seqtext = fread(fid,'*char')';
        fclose(fid);
    catch allExceptions %#ok<NASGU>
        error(message('bioinfo:seqread:CannotReadInput', filename));
    end
    [dispath, name, extension] = fileparts(filename); %#ok
    dispname = [name extension];
else
   error(message('bioinfo:seqread:CannotReadInput', filename));
end

% Guess the format from information in the first line
% Currently supported formats are 'FASTA','GENBANK','GENPEPT','EMBL','PDB'
% & raw text


firstLine = strread(seqtext,'%s','delimiter','\n');
firstLine = char(firstLine(find(~cellfun('isempty',firstLine),1))); %eliminate blank header lines

lw = lastwarn;

%  FASTA
if firstLine(1) == '>'
    try
        data = fastaread(seqtext);
        lastwarn(lw);
        filetype = 'FASTA';
        return
    catch allExceptions %#ok<NASGU>
    end
end

% GenBank or GenPept
if strncmpi('LOCUS',firstLine,5) == 1
    seqtext = char(strread(seqtext,'%s','delimiter','\n','whitespace',''));
    if ~isempty(strfind(firstLine,' aa '))
        try
            data = genpeptread(seqtext,'preambletext',getPreambleText);
            lastwarn(lw);
            filetype = 'GENPEPT';
            if isfield(data,'Sequence')
                return
            end
        catch allExceptions %#ok<NASGU>
        end
    end
    try
        data = genbankread(seqtext,'preambletext',getPreambleText);
        lastwarn(lw);
        filetype = 'GENBANK';
        if isfield(data,'Sequence')
            return
        end
    catch allExceptions %#ok<NASGU>
    end
end
% EMBL
if strncmpi('ID',firstLine,2) == 1

    try
        data = emblread(seqtext,'preambletext',getPreambleText);
        lastwarn(lw);
        filetype = 'EMBL';
        if isfield(data,'Sequence')
            return
        end
    catch allExceptions %#ok<NASGU>
    end

end
% PDB
if strncmpi('HEADER',firstLine,6) == 1
    try
        data = pdbread(seqtext);
        lastwarn(lw);
        filetype = 'HEADER';
        if isfield(data,'Sequence')
            return
        end
    catch allExceptions %#ok<NASGU>
    end
end

% Guess that the sequence is a single sequence in a raw text file
if isnt(firstLine) || isaa(firstLine)
    data.Header = sprintf('Sequence from %s',dispname);
    data.Sequence = regexprep(seqtext,'\W','');
    filetype = 'TEXT';
    return
end

error(message('bioinfo:seqread:UnknownFormat', dispname));

