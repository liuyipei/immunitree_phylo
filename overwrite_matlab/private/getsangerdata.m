function out=getsangerdata(key,varargin)
%GETSANGERDATA Retrieves sequence information from the PFAM database at the Sanger Institute.
%   OUT = GETSANGERDATA(KEY) queries for data identified with KEY in the
%   Sanger Institute databases, and returns a structure containing
%   information  for the protein family. Key may be a number or a string
%   containing the PFAM accession code.
%
%   GETSANGERDATA(...,'DATABASE',DB) will search in the specified 
%   database, DB, for data.  The accepted values are: 'align', 'hmm' or
%   'tree'. Default is 'align'.
%
%   GETSANGERDATA(...,'TYPE',TYPE) returns the alignments or phylogenetic
%   trees accessed by KEY. If TYPE='seed' only sequences used to generate
%   the HMM model, and if TYPE='full' all the sequences that hit the model
%   are returned. Valid only when database is 'align' or 'tree'. Default is
%   'full'.
%
%   GETSANGERDATA(...,'MODE',MODE) returns the Hidden Markov Model
%   identified by KEY. Global alignment model when MODE='ls' and local
%   alignment model when MODE='fs'. Valid only when database is 'hmm'.
%   Default is 'ls'.
%
%   GETSANGERDATA(...,'IGNOREGAPS',true) removes any gap symbol ('-' or
%   '.') from the sequences. Valid only when database is 'align'. Default
%   is false.  
%
%   GETSANGERDATA(...,'TOFILE',FILENAME) saves the data returned from the
%   database in the file FILENAME.
%
%   Examples:
%      
%   To retrieve an hmm profile model for global alignment to the 7-transmembrane 
%   receptor, Secretin family. (pfam key = PF00002)
%      hmmmodel  = getsangerdata(2,'database','hmm','mode','ls')
%   or
%      hmmmodel  = getsangerdata('PF00002','database','hmm','mode','ls')
%
%   To retrieve the aligned sequences, which generated such model.
%      pfamalign = getsangerdata(2,'database','align','type','seed')
%
%   To retrieve all other proteins aligned to this family. 
%      pfamalign = getsangerdata(2,'database','align','type','full')
%
%   For more information: http://pfam.sanger.ac.uk/
%
%   See also GETHMMALIGNMENT, GETHMMPROF, GETHMMTREE

% Copyright 2003-2005 The MathWorks, Inc.
% $Revision: 1.7.4.12 $   $Date: 2010/12/22 16:19:22 $

if ~usejava('jvm')
    error(message('bioinfo:getsangerdata:NeedJVM', mfilename));
end

% defaults
database = 'align';
type     = 'full' ;
mode     = 'ls';
mirrorSite = 'pfam.sanger.ac.uk';
tofile   = false;
ignoreGaps = false;
idType = 'acc';

if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:getsangerdata:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'database','type','mode','tofile','ignoregaps','mirror'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:getsangerdata:UnknownParameter', pname));
        elseif length(k)>1
            error(message('bioinfo:getsangerdata:AmbiguousParameter', pname)) ;
        else
            switch(k)
                 case 1  % database
                     okdbs = {'align','hmm','tree'} ;
                     val = find(strncmpi(pval,okdbs,length(pval)));
                     if isempty(val)
                         error(message('bioinfo:getsangerdata:InvalidDatabase'))
                     else
                         database = okdbs{val};
                     end
                 case 2 % type
                     okdbs = {'full','seed'} ;
                     val = find(strncmpi(pval,okdbs,length(pval)));
                     if isempty(val)
                         error(message('bioinfo:getsangerdata:InvalidType'))
                     else
                         type = okdbs{val} ;
                     end
                 case 3 % mode
                     okdbs = {'ls','fs'};
                     val = find(strncmpi(pval,okdbs,length(pval)));
                     if isempty(val)
                         error(message('bioinfo:getsangerdata:InvalidMode'))
                     else
                         mode = okdbs{val} ;
                     end
                 case 4  % tofile
                     if ischar(pval)
                         tofile = true;
                         filename = pval;
                     end
                 case 5  % ignore gaps
                    ignoreGaps = opttf(pval);
                    if isempty(ignoreGaps)
                        error(message('bioinfo:getsangerdata:IgnoreGapsNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 6
                    okdbs = {'sanger','janelia'};
                    val = find(strncmpi(pval,okdbs,length(pval)));
                    if isempty(val)
                        error(message('bioinfo:getsangerdata:InvalidMirrorSite'))
                    else
                        switch(val)
                            case 1
                                mirrorSite = 'pfam.sanger.ac.uk';
                            case 2
                                mirrorSite = 'pfam.janelia.org';
                        end
                    end
            end
        end
    end
end
% convert key to a string if it is a number
if isnumeric(key)
    str = num2str(key) ;
    key = 'PF00000' ;
    key(8-length(str):end)= str;
%Using HMM id for finding the model instead of accession number
elseif isempty(regexp(key,'PF\d{5,}','once'))
    idType = 'id';
end
% error if key isn't a string
if ~ischar(key)
    error(message('bioinfo:getsangerdata:NotString'))
end

% create the url that is used
switch database
    case 'align'
        searchurl = ['http://',mirrorSite,'/family/alignment/download/format?',idType,'=',key,'&format=fasta&alnType=',type];
    case 'hmm'
        searchurl = ['http://',mirrorSite,'/family/hmm?',idType,'=',key,'&mode=',mode];
    case 'tree'
        searchurl = ['http://',mirrorSite,'/family/tree/download?alnType=',type,'&',idType,'=',key];
end
% get the html file that is returned as a string
try
   s = urlread(searchurl) ;
catch allExceptions
   error(message('bioinfo:getsangerdata:ErrorDownloadingURL', searchurl))
end
% replace the html version of &
s=strrep(s,'&amp;','&');
% remove the start/end html markers
s=strrep(s,'<pre>','');
s=strrep(s,'</pre>', ' ' ) ;
s=strrep(s,'<PRE>','');
s=strrep(s,'</PRE>', ' ' ) ;
% pass data to respective parser, to create structure
switch database
    case 'align'
        if ignoreGaps
            out = fastaread(s,'ignoregaps',true);
        else
            out = fastaread(s);
        end
    case 'hmm'
        out = pfamhmmread(s);
    case 'tree'
        out = phytreeread(s);
end

%  write out file
if tofile == true
    fid = fopen(filename,'wt') ;
    if fid == (-1)
        error(message('bioinfo:getsangerdata:CouldNotOpenFile', filename));
    else
        fprintf(fid,'%s',s) ;
        fclose(fid);
    end
end
