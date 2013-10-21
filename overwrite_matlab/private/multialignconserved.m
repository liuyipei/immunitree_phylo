function conservedString = multialignconserved(alignment,varargin)
%MULTIALIGNCONSERVED creates a conservation string from a multiple alignment.
%
%   MULTIALIGNCONSERVED(ALIGNMENT) creates a string of the same length as
%   the sequenes in the alignment where '*' indicates positions with a
%   single, fully conserved residue, ':' indicates positions with strongly
%   conserved residues, and '.' indicates positions with weakly conserved
%   residues.
%
%   MULTIALIGNCONSERVED(..., 'STRONG',STRONGGROUPS) allows you to pass a cell
%   array of 'strong' groups. The default the 'strong' groups are:
%       {'STA', 'NEQK', 'NHQK', 'NDEQ', 'QHRK', 'MILV', 'MILF', 'HY',
%       'FYW'}
%
%   MULTIALIGNCONSERVED(..., 'WEAK',WEAKGROUPS) allows you to pass a cell
%   array of 'weak' groups. The default 'weak' groups are:
%   	{'CSA', 'ATV', 'SAG', 'STNK', 'STPA', 'SGND', 'SNDEQK', 'NDEQHK',
%   	'NEQHRK', 'FVLIM', 'HFY'}
%
%   Example:
%
%       % Reads the a multiple alignment of the gag polyprotein of several
%       % HIV strains.
%       gagaa = multialignread('aagag.aln')
%       multialignconserved(gagaa)
%
%   See also FASTAREAD, GETHMMALIGNMENT, MULTIALIGN, MULTIALIGNREAD,
%   MULTIALIGNVIEWER, SEQCONSENSUS, SEQDISP, SEQPROFILE.


%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2010/12/22 16:19:24 $

% Check inputs
bioinfochecknargin(nargin,1,mfilename);
[strongGroups,weakGroups] = parse_inputs(varargin{:});

% Convert the alignment sequence to an array that we can use to index
seqMatrix = double(aa2int(char({alignment.Sequence})));

% Ignore columns with gaps
gapCols = any(seqMatrix>20);

% Find perfect matches not counting any columns of all gaps
perfectCols = all(diff(seqMatrix) == 0,1);
perfectCols(perfectCols&gapCols) = false;

% Turn the cell arrays into matrices that we will index into
strongMatrix = cell2mat(cellfun(@makeIndexMatrix,strongGroups,'UniformOutput',false));
strongCols = false(size(gapCols));
strongOffsets = 20*(0:numel(strongGroups)-1);

% For columns that aren't perfect or gapped, see if they are strongly
% conserved.
unknownCols = ~(gapCols|perfectCols);
for i=find(unknownCols)
    strongCols(i) = any(all(strongMatrix(bsxfun(@plus,double(seqMatrix(:,i)),strongOffsets))));
end

% For columns that aren't perfect, gapped, or strongly conserved, see if they are weakly
% conserved.
weakMatrix = cell2mat(cellfun(@makeIndexMatrix,weakGroups,'UniformOutput',false));
weakCols = strongCols;
weakOffsets = 20*(0:numel(weakGroups)-1);
unknownCols = unknownCols&(~strongCols);
for i=find(unknownCols)
    weakCols(i) = any(all(weakMatrix(bsxfun(@plus,double(seqMatrix(:,i)),weakOffsets))));
end

% Now create the string and set the symbols for the hits.
conservedString = blanks(numel(gapCols));
conservedString(perfectCols) = '*';
conservedString(weakCols) = '.';
conservedString(strongCols) = ':';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function col = makeIndexMatrix(input)
col = zeros(20,1);
ndx = aa2int(input);
col(ndx) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [strongGroups,weakGroups] = parse_inputs(varargin)
% Parse the varargin parameter/value inputs

% Check that we have the right number of inputs
if rem(nargin,2) == 1
    error(message('bioinfo:multialignconserved:IncorrectNumberOfArguments', mfilename));
end

% The allowed inputs
okargs = {'strong','weak'};

% Set default values
strongGroups = {'STA',...
    'NEQK',...
    'NHQK',...
    'NDEQ',...
    'QHRK',...
    'MILV',...
    'MILF',...
    'HY',...
    'FYW'};

weakGroups = {'CSA',...
    'ATV',...
    'SAG',...
    'STNK',...
    'STPA',...
    'SGND',...
    'SNDEQK',...
    'NDEQHK',...
    'NEQHRK',...
    'FVLIM',...
    'HFY'};
% Loop over the values
for j=1:2:nargin
    % Lookup the pair
    [k, pval] = pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
    switch(k)
        case 1  % strong
            if ~isempty(pval)
                if ~iscellstr(pval)
                    error(message('bioinfo:multialignconserved:ConservedGroupNotCell'));
                end

                if ~all(cellfun(@isaa,pval))
                    error(message('bioinfo:multialignconserved:ConservedGroupNotCell'));
                end
                strongGroups = pval;
            end
        case 2  % weak
            if ~isempty(pval)
                if ~iscellstr(pval)
                    error(message('bioinfo:multialignconserved:SemiConservedGroupNotCell'));
                end

                if ~all(cellfun(@isaa,pval))
                    error(message('bioinfo:multialignconserved:ConservedGroupNotCell'));
                end
                weakGroups = pval;
            end
    end
end
