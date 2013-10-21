function out=getjaneliapfamdata(key,varargin)
%GETJANELIAPFAMDATA Retrieves sequence information from PFAM at Janelia Farms.
%   OUT = GETJANELIAPFAMDATA(KEY) queries for data identified with KEY in
%   the Janelia Farms PFAM databases, and returns a structure containing
%   information for the protein family. KEY may be a number or a string
%   containing the PFAM accession code.
%
%   GETJANELIAPFAMDATA(...,'DATABASE',DB) will search in the specified
%   database, DB, for data.  The accepted values are: 'align' or 'hmm'.
%   Default is 'align'.
%
%   GETJANELIAPFAMDATA(...,'TYPE',TYPE) returns the alignments accessed by
%   KEY. If TYPE='seed' only sequences used to generate the HMM model, and
%   if TYPE='full' all the sequences that hit the model are returned. Valid
%   only when database is 'align'. Default is 'full'.
%
%   GETJANELIAPFAMDATA(...,'MODE',MODE) returns the Hidden Markov Model
%   identified by KEY. Global alignment model when MODE='ls' and Local
%   alignment model when MODE='fs'. Valid only when database is 'hmm'.
%   Default is 'ls'.
%
%   GETJANELIAPFAMDATA(...,'IGNOREGAPS',true) removes any gap symbol ('-'
%   or '.') from the sequences. Valid only when database is 'align'.
%   Default is false.  
%
%   GETJANELIAPFAMDATA(...,'TOFILE',FILENAME) saves the data returned from
%   the database in the file FILENAME.
%
%   Examples:
%
%   To retrieve an hmm profile model for global alignment to the 7-transmembrane
%   receptor, Secretin family. (pfam key = PF00002)
%      hmmmodel  = getjaneliapfamdata(2,'database','hmm','mode','ls')
%   or
%      hmmmodel  = getjaneliapfamdata('PF00002.15','database','hmm','mode','ls')
%
%   To retrieve the aligned sequences, which generated such model.
%      pfamalign = getjaneliapfamdata(2,'database','align','type','seed')
%
%   To retrieve all other proteins aligned to this family.
%      pfamalign = getjaneliapfamdata(2,'database','align','type','full')
%
%   For more information:  http://pfam.janelia.org
%
%   See also GETHMMALIGNMENT, GETHMMPROF, GETHMMTREE

% Copyright 2003-2007 The MathWorks, Inc.
% $Revision: 1.1.6.6 $   $Date: 2010/12/22 16:19:20 $

% Stable links in the Janelia Farms PFAM site are described in
% http://pfam.janelia.org/help/stable_link.shtml 

%%%%% THIS FUNCTION IS RE-DIRECTED TO GETSANGERDATA. %%%%%%%%
% The reason for the redirect is that the two sites now run the exact same
% web-site, so script calling is exactly the same except for the base URL.
% This is handled in GETSANGERDATA.

out = getsangerdata(key,'mirror','janelia',varargin{:});

%%%%%%OLD CODE%%%%%%
% if ~usejava('jvm')
%     error('bioinfo:NeedJVM','%s requires Java.',mfilename);
% end
% 
% % defaults
% database = 'align';
% type     = 'full' ;
% mode     = 'ls';
% tofile   = false;
% ignoreGaps = false;
% idType = 'acc';
% 
% if nargin > 1
%     if rem(nargin,2) == 0
%         error('bioinfo:getjaneliapfamdata:IncorrectNumberOfArguments',...
%             'Incorrect number of arguments to %s.',mfilename);
%     end
%     okargs = {'database','type','mode','tofile','ignoregaps','mirror'};
%     for j=1:2:nargin-2
%         pname = varargin{j};
%         pval = varargin{j+1};
%         k = find(strncmpi(pname,okargs,numel(pname)));
%         if isempty(k)
%             error('bioinfo:getjaneliapfamdata:UnknownParameter',...
% 			'Unknown parameter name: %s.',pname);
%         elseif length(k)>1
%             error('bioinfo:getjaneliapfamdata:AmbiguousParameter',...
% 			'Ambiguous parameter name: %s.',pname) ;
%         else
%             switch(k)
%                 case 1  % database
%                     okdbs = {'align','hmm','tree'} ;
%                     val = find(strncmpi(pval,okdbs,length(pval)));
%                     if isempty(val)
%                         error('bioinfo:getjaneliapfamdata:InvalidDatabase','Invalid database name.')
%                     else
%                         database = okdbs{val};
%                     end
%                 case 2 % type
%                     okdbs = {'full','seed'} ;
%                     val = find(strncmpi(pval,okdbs,length(pval)));
%                     if isempty(val)
%                         error('bioinfo:getjaneliapfamdata:InvalidType','Invalid type.')
%                     else
%                         type = okdbs{val} ;
%                      end
%                  case 3 % mode
%                      okdbs = {'ls','fs'};
%                      val = find(strncmpi(pval,okdbs,length(pval)));
%                      if isempty(val)
%                          error('bioinfo:getjaneliapfamdata:InvalidMode','Invalid mode.')
%                      else
%                          mode = okdbs{val} ;
%                      end
%                  case 4  % tofile
%                      if ischar(pval)
%                          tofile = true;
%                          filename = pval;
%                      end
%                  case 5  % ignore gaps
%                     ignoreGaps = opttf(pval);
%                     if isempty(ignoreGaps)
%                         error('bioinfo:getjaneliapfamdata:IgnoreGapsNotLogical',...
%                               '%s must be a logical value, true or false.',...
%                               upper(char(okargs(k))));
%                     end
%                  case 6 % mirror is a ok arg, but it was already used in the 
%                         % wrapper function                                      
%             end
%         end
%     end
% end
% % convert key to a string if it is a number
% if isnumeric(key)
%     str = num2str(key) ;
%     key = 'PF00000' ;
%     key(8-length(str):end)= str;
% %Using HMM id for finding the model instead of accession number
% elseif isempty(regexp(key,'PF\d{5,}','once'))
%     idType = 'id';
% end
% % error if key isn't a string
% if ~ischar(key)
%     error('bioinfo:getjaneliapfamdata:NotString','Access Number is not a string.')
% end
% 
% % create the url that is used
% switch database
%     case 'align'
%         searchurl = ['http://pfam.janelia.org/family/alignment/download/format?',idType,'=',key,'&format=fasta&alnType=',type];
%     case 'hmm'
%         searchurl = ['http://pfam.janelia.org/family/gethmm?',idType,'=',key,'&mode=',mode];
%     case 'tree'
%         searchurl = ['http://pfam.janelia.org/family/tree/download?alnType=',type,'&',idType,'=',key];
% end
% % get the html file that is returned as a string
% try
%    s = urlread(searchurl) ;
% catch
%    error('bioinfo:getjaneliapfamdata:ErrorDownloadingURL',' Error downloading URL: %s', searchurl)
% end
% % replace the html version of &
% s=strrep(s,'&amp;','&');
% % remove the start/end html markers
% s=strrep(s,'<pre>','');
% s=strrep(s,'</pre>', ' ' ) ;
% s=strrep(s,'<PRE>','');
% s=strrep(s,'</PRE>', ' ' ) ;
% % pass data to respective parser, to create structure
% switch database
%     case 'align'
%         if ignoreGaps
%             out = fastaread(s,'ignoregaps',true);
%         else
%             out = fastaread(s);
%         end
%     case 'hmm'
%         out = pfamhmmread(s);
%     case 'tree'
%         out = phytreeread(s);
% end
% 
% %  write out file
% if tofile == true
%     fid = fopen(filename,'wt') ;
%     if fid == (-1)
%         error('bioinfo:getjaneliapfamdata:CouldNotOpenFile',...
%             'Unable to write to ''%s''.', filename);
%     else
%         fprintf(fid,'%s',s) ;
%         fclose(fid);
%     end
% end
