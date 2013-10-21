function gbout=getncbidata(accessnum,varargin)
% GETNCBIDATA retrieves sequence information from the NCBI databases.
%   GBOUT = GETNCBIDATA(ACCESSNUM)  searches for the Accession number in
%   the GenBank or GenPept database, and returns a structure containing
%   information for the sequence.
%
%   GBOUT = GETNCBIDATA(...,'DATABASE',DB) will search in the specified
%   database, DB, for data. The accepted values are: 'nucleotide' (default)
%   or 'protein'.
%
%   GBOUT = GETGENBANK(...,'PARTIALSEQ',SEQPARAMS) retrieves specified
%   subsequence for selected GenBank or GenPept file.  SEQPARAMS is a
%   two-element integer array containing the start and end positions of the
%   subsequence [START,END].  START is an integer between 1 and END and END
%   is an integer between START and length of sequence.
%
%   GBOUT = GETNCBIDATA(...,'TOFILE',FILENAME) saves the data returned from
%   the database in the file FILENAME using either the GenBank format or
%   the GenPept format.
%
%   GBOUT = GETNCBIDATA(...,'FILEFORMAT','FASTA') saves the data as a Fasta
%   formatted file. Additionally, the output structure GBOUT is a structure
%   with only the fields 'Header' and 'Sequence'.
%
%   GBOUT = GETNCBIDATA(...,'SEQUENCEONLY',true) returns just the sequence
%   as a character array. When the SEQUENCEONLY and TOFILE options are used
%   together, the output file is in the Fasta format.
%
%   Example:
%
%       nt = getncbidata('M10051','database','nucleotide')
%
%       aa = getncbidata('AAA59174','database','protein')
%
%   See: http://www.ncbi.nlm.nih.gov/About/disclaimer.html for more
%   information on NCBI data.
%
%   See also GENBANKREAD, GENPEPTREAD, GETGENBANK, GETGENPEPT.

% Copyright 2002-2011 The MathWorks, Inc.
% $Revision: 1.15.4.25 $   $Date: 2011/03/22 18:29:15 $

if ~usejava('jvm')
    error(message('bioinfo:getncbidata:NeedJVM', mfilename));
end
% default to 'nucleotide'
db = 'nucleotide';
dbfrm = '';
tofile = false;
seqonly = false;
fileformat = 'GenBank';
passThroughArgs = {};
partialseq = [];

if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:getncbidata:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'database','tofile','sequenceonly','fileformat',...
        'partialseq','preambletext','allfeatures'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:getncbidata:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:getncbidata:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % database
                    okdbs = {'nucleotide','protein'};
                    val = find(strncmpi(pval,okdbs,numel(pval)));
                    if length(val) == 1
                        dbfrs = {'gb','gp'};
                        db = okdbs{val};
                        dbfrm = dbfrs{val};
                    else
                        if isempty(val)
                            error(message('bioinfo:getncbidata:baddb'))
                        else
                            error(message('bioinfo:getncbidata:ambiguousdb', pval));
                        end
                    end
                case 2  % tofile
                    if ischar(pval)
                        tofile = true;
                        filename = pval;
                    end
                case 3  % sequenceonly
                    seqonly = opttf(pval);
                    if isempty(seqonly)
                        error(message('bioinfo:getncbidata:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 4 % fileformat
                    okformats = {'GenBank','GenPept','FASTA'};
                    val = find(strncmpi(pval,okformats,numel(pval)));
                    if length(val) == 1
                        fileformat = okformats{val};
                        dbfrm = '';
                    else
                        if isempty(val)
                            error(message('bioinfo:getncbidata:badformat'));
                        else
                            error(message('bioinfo:getncbidata:ambiguousfileformat', pval));
                        end
                    end
                case 5 %partialseq
                    if iscell(pval) || ischar(pval)
                        error('bioinfo:getncbidata:PartialseqWrongTypes',...
                            ['PARTIALSEQ parameters must be a two-element'...
                            ' integer array, not a cell or character array']);
                    elseif numel(pval)==1
                        error(message('bioinfo:getncbidata:TooFewPartialseq'));
                    elseif numel(pval)>2
                        error(message('bioinfo:getncbidata:TooManyPartialseq'));
                    else
                        partialseq = pval;
                    end
                    
                case {6,7}  % 'preambletext'
                    passThroughArgs{end+1} = pname; %#ok<AGROW>
                    passThroughArgs{end+1} = opttf(pval); %#ok<AGROW>
            end
        end
    end
end

% if a fileformat was specified as an input, then dbfrm should now be empty.  Use the
% specified format
if isempty(dbfrm)
    if strcmpi(fileformat,'GenBank')
        dbfrm = 'gb';
    elseif strcmpi(fileformat,'GenPept')
        dbfrm = 'gp';
    elseif strcmpi(fileformat,'Fasta')
        dbfrm = 'fasta';
    end
end

% if only the sequence is required, bring the data as a FASTA formatted
% file
if seqonly
    fileformat = 'Fasta';
    dbfrm = 'fasta';
end

% convert accessnum to a string if it is a number
if isnumeric(accessnum)
    accessnum = num2str(accessnum);
end

% error if accessnum isn't a string
if ~ischar(accessnum)
    error(message('bioinfo:getncbidata:NotString'))
end

% error if START and END of partialseq are not ordered.
if ~isempty(partialseq) && partialseq(2)<partialseq(1)
    error(message('bioinfo:getncbidata:EndLessThanStart', partialseq( 1 ), partialseq( 2 )));
end


%PDB entries in Protein have an underscore in the accession number in the
%record, but the accession number to retrieve the file cannot have the
%underscore or no file will be found. (eg. 1JTS_O in file, but query with
%1JTSO).  This will remove the underscore from the accession number.
%Accession number must have the format 4 characters_1 character.
if ~isempty(regexp(accessnum,'\<\w{4}_\w{1}\>','once'))&&strcmpi(db,'protein')
    accessnum = strrep(accessnum,'_','');
end


%Subfunction ACCESSION2GI does the querying and searching of the accession
%number in specified databases to get the GI numbers.
%
%The first time ACCESSION2GI is called it performs the query and search
%with the tag [Accession] to expedite the query/search.
[giID,db] = accession2gi(accessnum,db,'quick');

%The second time ACCESSION2GI is called it performs the query and search
%WITHOUT the tag [Accession].  This is slower, but finds files that are
%dead, replaced or suppressed so an appropriate warning can be thrown and
%the file retrieved.
if isempty(giID)
    [giID,db] = accession2gi(accessnum,db,'slow');
    if isempty(giID)
        error('bioinfo:getncbidata:NoItemsFound',...
            ['The key %s was not found in the %s database at this time.\n',...
            'Please check that the input is a valid accession number or try again.\n\n',...
            'NOTE:  This function is dependent on NCBI''s Entrez tools and sequence databases. ',...
            'Changes to either may cause this function to break.  Please check Bug Report 492836 at ',...
            'http://www.mathworks.com/support/bugreports/details.html?rp=492836 for any fixes or ',...
            'updates to this function due to changes in NCBI''s Entrez tools or sequence databases.'],...
            accessnum,db) ;
    end
end

% create url to retrieve record as text
if isempty(partialseq)
    retrieveurls = ...
        ['http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=' char(db)...
        '&id=' char(giID{1}) '&rettype=',dbfrm,'&retmode=text'];
else
    if numel(partialseq)==2
        retrieveurls = ...
            ['http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=' char(db)...
            '&id=' char(giID{1}) '&rettype=',dbfrm,'&retmode=text'...
            '&seq_start=',num2str(partialseq(1)),'&seq_stop=',num2str(partialseq(2))];
    else
        retrieveurls = ...
            ['http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=' char(db)...
            '&id=' char(giID{1}) '&rettype=',dbfrm,'&retmode=text'...
            '&seq_start=',num2str(partialseq(1)),'&seq_stop=',num2str(partialseq(2)),...
            '&strand=',num2str(partialseq(3))];
    end
end

% shortcut to display hyperlinks without actually downloading the data
if (tofile == false) && (nargout==0)
    SearchURL = ['http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=' db '&id=' accessnum];
    disp([char(9) 'SearchURL: <a href="' SearchURL '">' accessnum '</a> -> ' SearchURL]);
    disp([char(9) 'RetrieveURL: <a href="' retrieveurls '">' char(giID{1}) '</a> -> ' retrieveurls]);
    return
end

% determine the file where URLWRITE will download the data
if (tofile == true) &&  ~exist(filename,'file')
    dataIsInTemp = false;
    urlwritefilename = filename;
else
    dataIsInTemp = true;
    urlwritefilename = tempname;
end

[~,status] = urlwrite(retrieveurls,urlwritefilename);
if ~status
    error(message('bioinfo:getncbidata:ErrorURLWrite', urlwritefilename));
end

% check if downloaded file is empty
filestatus = dir(urlwritefilename);
if  filestatus.bytes == 0
    error(message('bioinfo:getncbidata:ErrorFetchingData', accessnum, db));
end

% check if sequence no longer exists in the database
fid = fopen(urlwritefilename,'rt');
if fid == (-1)
    error('bioinfo:getncbidata:CouldNotOpenFile',...
        'Unable to read to ''%s''.', urlwritefilename);
end
genbankdata = fscanf(fid,'%c',16384);
fclose(fid);
if ~isempty(regexp(genbankdata,'Error: download dataset is empty','match'))
    error('bioinfo:getncbidata:RecordDiscontinued',...
        'File %s has been replace or removed from %s database',accessnum,db);
end

if nargout
    % read new file with GENBANKREAD or FASTAREAD to create structure
    try
        switch dbfrm
            case 'gb'
                gbout = genbankread(urlwritefilename,passThroughArgs{:});
            case 'gp'
                gbout = genpeptread(urlwritefilename,passThroughArgs{:});
            case 'fasta'
                gbout = fastaread(urlwritefilename);
        end
    catch allExceptions
        error(message('bioinfo:getncbidata:URLreturnedIncorrectData', filename))
        % this error is random, depends on the reliability of the Internet connection
    end
    
    % Append URLs to output structure
    if ~strcmpi(fileformat,'FASTA')
        gbout.SearchURL = ['http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=' db '&id=' accessnum];
        gbout.RetrieveURL = retrieveurls;
    end
    
    if seqonly == true
        gbout = gbout.Sequence;
    end
end

% append urlwritefilename to filename when necessary
if (tofile == true) && dataIsInTemp
    % The default behavior is to append to the file. If the file exists we
    % warn, however, we don't want to warn about appending until we are
    % sure we can open the file.
    fid1 = fopen(filename,'at');
    if fid1 == (-1)
        error('bioinfo:getncbidata:CouldNotOpenFile',...
            'Unable to write to ''%s''.', filename);
    end
    c1 = onCleanup(@()fclose(fid1));
    warning(message('bioinfo:getncbidata:AppendToFile', filename));
    fid2 = fopen(urlwritefilename,'rt');
    if fid2 == (-1)
        error('bioinfo:getncbidata:CouldNotOpenFile',...
            'Unable to read to ''%s''.', urlwritefilename);
    end
    c2 = onCleanup(@()fclose(fid2));
    blocksize = 16777216;
    count = 1;
    % copy data by blocks
    while count
        [genbankdata count] = fscanf(fid2,'%c',blocksize);
        fprintf(fid1,'%c',genbankdata);
    end
end

% delete temporal file used by urlwrite
if dataIsInTemp
    clear c2; % closes file
    delete(urlwritefilename)
end

function [giID,db] = accession2gi(accessnum,db,qsFlag)
%This subfunction preforms an esearch and esummary to find given
%accession numbers in the protein or nucleotide databases.  It returns the
%giID for the accession number.  qsFlag can be 'quick' or 'slow'.

% In the abcense of GQUERY, we try sequentially 4 databases.
if strcmpi('nucleotide',db)
    % Try all four databases starting with the nucleotide databases.
    trydb = {'nuccore','nucgss','nucest','protein'};
    typedb = {'nucleotide','nucleotide','nucleotide','protein'};
else
    % Try all four databases starting with the protein databases.
    trydb = {'protein','nuccore','nucgss','nucest',};
    typedb = {'protein','nucleotide','nucleotide','nucleotide'};
end

for i = 1:4
    
    db_i = trydb{i};
    type_i = typedb{i};
    
    %Build esearch URL
    if strcmpi(qsFlag,'quick')
        %Using the [Accession] key word expedites the query
        searchurl = ...
            ['http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db='...
            db_i '&term=' accessnum '[Accession]'];
    else
        %Not using the [Accession] key word slows the query speed, but allows
        %dead, suppressed, and replaced files to be found.
        searchurl = ...
            ['http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db='...
            db_i '&term=' accessnum];
    end
    
    %use esearch to get GI ID for specified accession number
    searchXML = urlread(searchurl);
    
    if isempty(searchXML)
        error(message('bioinfo:getncbidata:EsearchProblem', db_i));
    end
    
    giID = regexp(searchXML,'<Id>(\w+)</Id>','tokens');
    
    if ~isempty(giID)
        break
    end
end



%error if accession number not found in specified database
if isempty(giID)
    %do not error, return empty string to be able to toggle between fast
    %and slow query/search, the main function will error if an empty string
    %is returned in both tries
    giID = '';
    return;
end

%Error if sequence was found in other database
if ~strcmpi(db,type_i)
    if strcmpi(type_i,'protein')
        program = 'getgenpept';
    else
        program = 'getgenbank';
    end
    error('bioinfo:getncbidata:IncorrectKnownDatabaseSuggest',...
        ['The key %s was not found in the %s database,\n'...
        'but it is in the %s database.  Try using: %s(''%s'')'],...
        accessnum,db,'Nucleotide',program,accessnum);
end

db = db_i;

%use esummary to get DocSum for record and parse XML to get record
%status, replacement record and comments
summaryurl = ['http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=' db '&id=' char(giID{1})];
summaryXML = urlread(summaryurl);

if isempty(summaryXML)
    error(message('bioinfo:getncbidata:EsummaryProblem', accessnum, db));
end

sumExp = ['<Item Name="Caption" Type="String">(?<accnum>\w*)</Item>.*?'...
    '<Item Name="Status" Type="String">(?<status>\w*)</Item>\s*'...
    '<Item Name="ReplacedBy" Type="String">(?<replaced>.*?)</Item>\s*'...
    '<Item Name="Comment" Type="String"><!\[CDATA\[\s*(?<comment>(.*?))\]\]></Item>'];

recordStatus = regexp(summaryXML,sumExp,'names');

%Determine if a record is dead, suppressed or replaced and warn
%appropriately
if isempty(recordStatus)
    %Error if nothing was parsed out of the summary report
    error(message('bioinfo:getncbidata:EsummaryResultsProblem', accessnum, db));
elseif strcmpi(recordStatus.status,'dead')
    %Warn if record has been discontinued
    warning('bioinfo:getncbidata:RecordDiscontinued',...
        'The record %s has been discontinued in the NCBI database.',...
        accessnum);
elseif strcmpi(recordStatus.status,'suppressed')
    %Warn if record has been suppressed
    warning('bioinfo:getncbidata:RecordSuppressed',...
        ['The record %s has been suppressed in the NCBI database.\n'...
        'Returning record %s'],...
        accessnum,accessnum);
    if ~isempty(recordStatus.replaced)
        disp(['NOTE: Record ',accessnum,' has been replaced by ',recordStatus.replaced]);
    end
elseif strcmpi(recordStatus.status,'replaced')
    %Warn if record has been replaced
    warning('bioinfo:getncbidata:RecordReplaced',...
        ['The record %s has been replaced by record %s\nin the %s database.',...
        '  Returning record %s.'],...
        accessnum,recordStatus.replaced,db,accessnum);
elseif ~strcmpi(recordStatus.accnum,accessnum)
    %Warn if record has been replaced and new record is automatically
    %returned (eg. getgenbank('M73485'))
    warning('bioinfo:getncbidata:RecordAutoReplaced',...
        ['The record %s has been replaced by %s.\n',...
        'Returning record %s'],...
        accessnum,recordStatus.accnum,char(giID{1}));
end


