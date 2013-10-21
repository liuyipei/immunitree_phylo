function [score, alignment, startat, matrices] = nwalign(seq1,seq2,varargin)
%NWALIGN performs Needleman-Wunsch global alignment of two sequences.
%
%   NWALIGN(SEQ1, SEQ2) returns the score (in bits) for the optimal
%   alignment. Note: The scale factor used to calculate the score is
%   provided by the scoring matrix info (see below). If this is not
%   defined, then NWALIGN returns the raw score.
%
%   [SCORE, ALIGNMENT] = NWALIGN(SEQ1, SEQ2) returns a string showing an
%   optimal global alignment of amino acid (or nucleotide) sequences SEQ1
%   and SEQ2.
%
%   [SCORE, ALIGNMENT, STARTAT] = NWALIGN(SEQ1, SEQ2)  returns a 2x1 vector
%   with the starting point indices indicating the starting point of the
%   alignment in the two sequences. Note: this output is for consistency
%   with SWALIGN and will always be [1;1] because this is a global
%   alignment.
%
%   NWALIGN(..., 'ALPHABET', A) specifies whether the sequences are
%   amino acids ('AA') or nucleotides ('NT'). The default is AA.
%
%   NWALIGN(..., 'SCORINGMATRIX', matrix) defines the scoring matrix to be
%   used for the alignment. The default is BLOSUM50 for AA or NUC44 for NT.
%
%   NWALIGN(..., 'SCALE' ,scale) indicates the scale factor of the scoring
%   matrix to return the score using arbitrary units. If the scoring matrix
%   Info also provides a scale factor, then both are used.
%
%   NWALIGN(..., 'GAPOPEN', penalty) defines the penalty for opening a gap
%   in the alignment. The default gap open penalty is 8.
%
%   NWALIGN(..., 'EXTENDGAP', penalty) defines the penalty for extending a
%   gap in the alignment. If EXTENDGAP is not specified, then extensions to
%   gaps are scored with the same value as GAPOPEN.
%
%   NWALIGN(..., 'SHOWSCORE', true) displays the scoring space and the
%   winning path.
%
%
%   Examples:
%
%       % Return the score in bits and the global alignment using the
%       % default scoring matrix (BLOSUM50).
%       [score, align] = nwalign('VSPAGMASGYD', 'IPGKASYD')
%
%       % Use user-specified scoring matrix and "gap open" penalty.
%       [score, align] = nwalign('IGRHRYHIGG', 'SRYIGRG',...
%                               'scoringmatrix', @pam250, 'gapopen',5)
%
%       % Return the score in nat units (nats).
%       [score, align] = nwalign('HEAGAWGHEE', 'PAWHEAE', 'scale', log(2))
%
%       % Display the scoring space and the winning path.
%       nwalign('VSPAGMASGYD', 'IPGKASYD', 'showscore', true)
%
%   See also ALIGNDEMO, BLOSUM, LOCALALIGN, MULTIALIGN, NT2AA, PAM,
%   PROFALIGN, SEQDOTPLOT, SHOWALIGNMENT, SWALIGN.

%   References:
%   R. Durbin, S. Eddy, A. Krogh, and G. Mitchison. Biological Sequence
%   Analysis. Cambridge UP, 1998.
%   Needleman, S. B., Wunsch, C. D., J. Mol. Biol. (1970) 48:443-453

%   Copyright 2002-2008 The MathWorks, Inc.
%   $Revision: 1.22.6.18 $  $Date: 2009/05/07 18:15:33 $


gapopen = -8;
gapextend = -8;
setGapExtend = false;
showscore=false;
isAminoAcid = true;
scale=1;

% If the input is a structure then extract the Sequence data.
if isstruct(seq1)
    seq1 = seqfromstruct(seq1);
end
if isstruct(seq2)
    seq2 = seqfromstruct(seq2);
end
if nargin > 2
    if rem(nargin,2) == 1
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'scoringmatrix','gapopen','extendgap','alphabet','scale','showscore'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % scoring matrix
                    if isnumeric(pval)
                        ScoringMatrix = pval;
                    else
                        if ischar(pval)
                            pval = lower(pval);
                        end
                        try
                            [ScoringMatrix,ScoringMatrixInfo] = feval(pval);
                        catch allExceptions
                            error('Bioinfo:InvalidScoringMatrix','Invalid scoring matrix.');
                        end
                    end
                case 2 %gap open penalty
                    gapopen = -pval;
                case 3 %gap extend penalty
                    gapextend = -pval;
                    setGapExtend = true;
                case 4 %if sequence is nucleotide
                    if strcmpi(pval,'nt')
                        isAminoAcid = false;
                    end
                case 5 % scale
                    scale=pval;
                case 6 % showscore
                    showscore = pval == true;
            end
        end
    end
end

% setting the default scoring matrix
if ~exist('ScoringMatrix','var')
    if isAminoAcid
        [ScoringMatrix,ScoringMatrixInfo] = blosum50;
    else
        [ScoringMatrix,ScoringMatrixInfo] = nuc44;
    end
end


% getting the scale from ScoringMatrixInfo, if it exists
if exist('ScoringMatrixInfo','var') && isfield(ScoringMatrixInfo,'Scale')
    scale=scale*ScoringMatrixInfo.Scale;
end

% handle properly "?" characters typically found in pdb files
if isAminoAcid
    if ischar(seq1)
        seq1 = strrep(seq1,'?','X');
    else
        seq1(seq1 == 26) = 23;
    end
    if ischar(seq2)
        seq2 = strrep(seq2,'?','X');
    else
        seq2(seq2 == 26) = 23;
    end
end

% check input sequences
% commented by Joni
% if isAminoAcid && ~(isaa(seq1) && isaa(seq2))
%     error('Bioinfo:InvalidAminoAcidSequences',...
%         'Both sequences must be amino acids, use ALPHABET = ''NT'' for aligning nucleotides.');
% elseif ~isAminoAcid && ~(isnt(seq1) && isnt(seq2))
%     error('Bioinfo:InvalidNucleotideSequences',...
%         'When ALPHABET = ''NT'', both sequences must be nucleotides.');
% end

% use numerical arrays for easy indexing
if ischar(seq1)
    seq1=upper(seq1); %the output alignment will be all uppercase
    if isAminoAcid
        intseq1 = aa2int(seq1);
    else
        intseq1 = nt2int(seq1);
    end
else
    intseq1 = uint8(seq1);
    seq1 = intseq1;
% commented by Joni, March 11th, 2010    
%    if isAminoAcid
%        seq1 = int2aa(intseq1);
%    else
%        seq1 = int2nt(intseq1);
%    end
end
if ischar(seq2)
    seq2 = upper(seq2); %the output alignment will be all uppercase
    if isAminoAcid
        intseq2 = aa2int(seq2);
    else
        intseq2 = nt2int(seq2);
    end
else
    intseq2 = uint8(seq2);
    seq2 = intseq2;
% commented by Joni, March 11th, 2010    
%     if isAminoAcid
%         seq2 = int2aa(intseq2);
%     else
%         seq2 = int2nt(intseq2);
%     end
end


m = length(seq1);
n = length(seq2);
if ~n||~m
    error('Bioinfo:InvalidLengthSequences','Length of input sequences must be greater than 0');
end

% If unknown, ambiguous or gaps appear in the sequence, we need to make
% sure that ScoringMatrix can handle them.

% possible values are
% B  Z  X  *  -  ?
% 21 22 23 24 25 26

scoringMatrixSize = size(ScoringMatrix,1);

highestVal = max([intseq1, intseq2]);
if highestVal > scoringMatrixSize
    % if the matrix contains the 'Any' we map to that
    if isAminoAcid
        anyVal = aa2int('X');
    else
        anyVal = nt2int('N');
    end
    if scoringMatrixSize >= anyVal
        intseq1(intseq1>scoringMatrixSize) = anyVal;
        intseq2(intseq2>scoringMatrixSize) = anyVal;
    else
        error('Bioinfo:InvalidSymbolsInInputSequences',...
            'Sequences contain symbols that cannot be handled by the given scoring matrix.');
    end
end

if setGapExtend 
    % flip order of input sequences for consistency with older versions
    if showscore % return the score matrices
        [score, path(:,2), path(:,1), F(:,:,1), F(:,:,2), F(:,:,3)] = ...
            affinegapmex(intseq2, intseq1, gapopen, gapextend, ScoringMatrix, 1); 
        
    elseif nargout ==  4 % return score matrices and pointer matrices
        [score, path(:,2), path(:,1), F(:,:,1), F(:,:,2), F(:,:,3), pointer] = ...
            affinegapmex(intseq2, intseq1, gapopen, gapextend, ScoringMatrix, 1); 
        pointer = shiftdim(pointer,1); % for backward compatibility
           
    else % return only score and alignment
        [score, path(:,2), path(:,1)] = affinegapmex(intseq2, intseq1, ...
            gapopen, gapextend, ScoringMatrix, 1); 
    end

else
     % flip order of input sequences for consistency with older versions
    if showscore % return the score matrices
        [score, path(:,2), path(:,1), F] = ...
            simplegapmex(intseq2, intseq1, gapopen, ScoringMatrix, 1);
        
    elseif nargout == 4
        [score, path(:,2), path(:,1), F, pointer] = ...
            simplegapmex(intseq2, intseq1, gapopen, ScoringMatrix, 1);
        
    else
        [score, path(:,2), path(:,1)] = ...
            simplegapmex(intseq2, intseq1, gapopen, ScoringMatrix, 1);
    end
end

path = path(sum(path,2)>0,:);
path = flipud(path);

% re-scaling the output score
score = scale * score;

if nargout<=1 && ~showscore
    return
end

% setting the size of the alignment
alignment = repmat(('- -')',1,size(path,1));

% adding sequence to alignment
alignment(1,path(:,1)>0) = '0'+intseq1;
alignment(3,path(:,2)>0) = '0'+intseq2;
% commented by Joni, March 11th, 2010    
% alignment(1,path(:,1)>0) = seq1;
% alignment(3,path(:,2)>0) = seq2;

% find locations where there are no gaps
h=find(all(path>0,2));
if isAminoAcid
    noGaps1=aa2int(alignment(1,h));
    noGaps2=aa2int(alignment(3,h));
else
    noGaps1=alignment(1,h);
    noGaps2=alignment(3,h);
% commented by Joni, March 11th, 2010    
%     noGaps1=nt2int(alignment(1,h));
%     noGaps2=nt2int(alignment(3,h));
end

% erasing symbols that cannot be scored
htodel=max([noGaps1;noGaps2])>scoringMatrixSize;
h(htodel)=[];
noGaps1(htodel)=[];
noGaps2(htodel)=[];

% score pairs with no gap
value = ScoringMatrix(sub2ind(size(ScoringMatrix),double(noGaps1),double(noGaps2)));

% insert symbols of the match string into the alignment
alignment(2,h(value>=0)) = ':';
alignment(2,h(noGaps1==noGaps2)) = '|';

startat = [1;1];

% undocumented fourth output -- score and pointer matrices
if nargout > 3
    matrices.Scores = F;
    matrices.Pointers = pointer;
end

if showscore
    figure
    F=scale.*max(F(2:end,2:end,:),[],3);
    clim=max(max(max(abs(F(~isinf(F))))),eps);
    imagesc(F,[-clim clim]);
    colormap(privateColorMap(1));
    set(colorbar,'YLim',[min([F(:);-eps]) max([F(:);eps])])
    title('Scoring Space and Winning Path')
    xlabel('Sequence 1')
    ylabel('Sequence 2')
    hold on
    plot(path(all(path>0,2),1),path(all(path>0,2),2),'k.')
end

%=== SIMPLEGAP is now a mex function ===%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [F, pointer] = simplegap(intseq1,m,intseq2,n,ScoringMatrix,gap)
% % Standard Needleman-Wunsch algorithm
% 
% % set up storage for dynamic programming matrix
% F = zeros(n+1,m+1);
% F(2:end,1) = gap * (1:n)';
% F(1,2:end) = gap * (1:m);
% 
% % and for the back tracing matrix
% pointer= repmat(uint8(4),n+1,m+1);
% pointer(:,1) = 2;  % up
% pointer(1,1) = 1;
% 
% 
% % initialize buffers to the first column
% ptr = pointer(:,2); % ptr(1) is always 4
% currentFColumn = F(:,1);
% 
% % main loop runs through the matrix looking for maximal scores
% for outer = 2:m+1
% 
%     % score current column
%     scoredMatchColumn = ScoringMatrix(intseq2,intseq1(outer-1));
%     % grab the data from the matrices and initialize some values
%     lastFColumn    = currentFColumn;
%     currentFColumn = F(:,outer);
%     best = currentFColumn(1);
% 
%     for inner = 2:n+1
%         % score the three options
%         up       = best + gap;
%         left     = lastFColumn(inner) + gap;
%         diagonal = lastFColumn(inner-1) + scoredMatchColumn(inner-1);
% 
%         % max could be used here but it is quicker to use if statements
%         if up > left
%             best = up;
%             pos = 2;
%         else
%             best = left;
%             pos = 4;
%         end
% 
%         if diagonal >= best
%             best = diagonal;
%             ptr(inner) = 1;
%         else
%             ptr(inner) = pos;
%         end
%         currentFColumn(inner) = best;
% 
%     end % inner
%     % put back updated columns
%     F(:,outer)   = currentFColumn;
%     % save columns of pointers
%     pointer(:,outer)  = ptr;
% end % outer
% 


%=== AFFINEGAP is now a mex function ===%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [F,pointer] = affinegap(intseq1,m,intseq2,n,ScoringMatrix,gapopen,gapextend)
% % Needleman-Wunsch algorithm modified to handle affine gaps
% 
% % Set states
% inAlign =   1;
% inGapUp =   2;
% inGapLeft = 3;
% numStates = 3;
% 
% % Set up storage for dynamic programming matrix:
% % for keeping the maximum scores for every state
% 
% F =  zeros(n+1,m+1,numStates);
% F(:,1,:) = -inf;
% F(1,:,:) = -inf;
% F(1,1,inAlign) = 0;
% 
% F(2:end,1,inGapUp)   = gapopen + gapextend * (0:n-1)';
% F(1,2:end,inGapLeft) = gapopen + gapextend * (0:m-1);
% 
% % and for the back tracing pointers
% pointer(n+1,m+1,numStates) = uint8(0);
% pointer(2:end,1,inGapUp)   = 2;  % up
% pointer(1,2:end,inGapLeft) = 4;  % left
% 
% % initialize buffers to the first column
% ptrA = pointer(:,1,inAlign);
% ptrU = pointer(:,1,inGapLeft);
% ptrL = pointer(:,1,inGapUp);
% 
% currentFColumnA = F(:,1,inAlign);
% currentFColumnU = F(:,1,inGapUp);
% currentFColumnL = F(:,1,inGapLeft);
% 
% % main loop runs through the matrix looking for maximal scores
% for outer = 2:m+1
%     % score current column
%     scoredMatchColumn = ScoringMatrix(intseq2,intseq1(outer-1));
%     % grab the data from the matrices and initialize some values for the
%     % first row the most orderly possible
%     lastFColumnA    = currentFColumnA;
%     currentFColumnA = F(:,outer,inAlign);
%     bestA           = currentFColumnA(1);
%     currentinA      = lastFColumnA(1);
% 
%     lastFColumnU    = currentFColumnU;
%     currentFColumnU = F(:,outer,inGapUp);
%     bestU           = currentFColumnU(1);
% 
%     lastFColumnL    = currentFColumnL;
%     currentFColumnL = F(:,outer,inGapLeft);
%     currentinGL     = lastFColumnL(1);
% 
%     for inner = 2:n+1
% 
%         % grab the data from the columns the most orderly possible
%         upOpen      = bestA + gapopen;
%         inA         = currentinA;
%         currentinA  = lastFColumnA(inner);
%         leftOpen    = currentinA + gapopen;
% 
%         inGL        = currentinGL;
%         currentinGL = lastFColumnL(inner);
%         leftExtend  = currentinGL + gapextend;
% 
%         upExtend = bestU + gapextend;
%         inGU     = lastFColumnU(inner-1);
% 
%         % operate state 'inGapUp'
% 
%         if upOpen > upExtend
%             bestU = upOpen; ptr = 1;   % diagonal
%         elseif upOpen < upExtend
%             bestU = upExtend; ptr = 2; % up
%         else % upOpen == upExtend
%             bestU = upOpen; ptr = 3;   % diagonal and up
%         end
%         currentFColumnU(inner)=bestU;
%         ptrU(inner)=ptr;
% 
%         % operate state 'inGapLeft'
% 
%         if leftOpen > leftExtend
%             bestL = leftOpen; ptr = 1;   % diagonal
%         elseif leftOpen < leftExtend
%             bestL = leftExtend; ptr = 4; % left
%         else % leftOpen == leftExtend
%             bestL = leftOpen; ptr = 5;   % diagonal and left
%         end
%         currentFColumnL(inner) = bestL;
%         ptrL(inner) = ptr;
% 
%         % operate state 'inAlign'
% 
%         if  inA > inGU
%             if inA > inGL
%                 bestA = inA; ptr = 1;  % diagonal
%             elseif inGL > inA
%                 bestA = inGL; ptr = 4; % left
%             else
%                 bestA = inA; ptr = 5;  % diagonal and left
%             end
%         elseif inGU > inA
%             if inGU > inGL
%                 bestA = inGU; ptr = 2; % up
%             elseif inGL > inGU
%                 bestA = inGL; ptr = 4; % left
%             else
%                 bestA = inGU; ptr = 6; % up & left
%             end
%         else
%             if inA > inGL
%                 bestA = inA; ptr = 3;  % diagonal & up
%             elseif inGL > inA
%                 bestA = inGL; ptr = 4; % left
%             else
%                 bestA = inA; ptr = 7;  % all
%             end
%         end
% 
%         bestA = bestA + scoredMatchColumn(inner-1);
%         currentFColumnA(inner) = bestA;
%         ptrA(inner) = ptr;
% 
%     end %inner
% 
%     % put back updated columns
%     F(:,outer,inGapLeft) = currentFColumnL;
%     F(:,outer,inGapUp)   = currentFColumnU;
%     F(:,outer,inAlign)   = currentFColumnA;
%     % save columns of pointers
%     pointer(:,outer,inAlign)  = ptrA;
%     pointer(:,outer,inGapUp)  = ptrU;
%     pointer(:,outer,inGapLeft)= ptrL;
% end %outer
end

function pcmap = privateColorMap(selection)
%PRIVATECOLORMAP returns a custom color map
switch selection
    case 1, pts = [0 0 .3 20;
            0 .1 .8 25;
            0 .9 .5 15;
            .9 1 .9 8;
            1 1 0 26;
            1 0 0 26;
            .4 0 0 0];
    otherwise, pts = [0 0 0 128; 1 1 1 0];
end
xcl=1;
for i=1:size(pts,1)-1
    xcl=[xcl,i+1/pts(i,4):1/pts(i,4):i+1]; %#ok<AGROW>
end
pcmap = interp1(pts(:,1:3),xcl);

end

function result = isaa(seq)
%ISAA True for amino acid sequences.
%   ISAA(SEQ) returns 1 for a amino acid sequence, 0 otherwise. Valid
%   symbols are A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,B,Z,X and *.  
%
%   See also ISDNA, ISRNA, ISNT.

%   Copyright 2002-2004 The MathWorks, Inc.
%   $Revision: 1.8.6.6 $  $Date: 2005/06/09 21:57:23 $

persistent maxval
if isempty(maxval)
    [dummy, map] = aa2int('a'); %#ok
    maxval = max(map);
end
if ischar(seq)
    try
        seq = aa2int(seq);
    catch
        result = false;
        return
    end
end

result = ~any(any(seq <= 0 | seq > maxval | seq ~= floor(seq)));

%exp = '[^ARNDCQEGHILKMFPSTWYVBZX]';
%result = isempty(regexpi(seq,exp,'once'));
end

function result = isnt(seq,varargin)
%ISNT True for nucleotide sequences.
%   ISNT(SEQ) returns 1 for a DNA sequence, 0 otherwise. Valid symbols are
%   A,C,G,T,U,N,R,Y,K,M,S,W,B,D,H,V and *.
%
%   ISNT(...,'ACGTUOnly',true) returns 1 only if the sequence contains
%   A,C,G and T or U only.   
%
%   See also ISDNA, ISRNA, ISAA.

%   Copyright 2002-2004 The MathWorks, Inc.
%   $Revision: 1.10.6.8 $  $Date: 2007/01/12 01:27:23 $

acgtonly = false;
if nargin > 1
    if rem(nargin,2)== 0
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'acgtuonly'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
			'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbiguousParameterName',...
			'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % others forces everything except ACTG to be the unknown value
                    if islogical(pval)
                        acgtonly = pval;
                    else
                        error('Bioinfo:InvalidACGTOnly',...
                            'Invalid value for ACGTonly parameter. Value must be true or false.');
                    end
            end
        end
    end
end

persistent maxval
    if isempty(maxval)
        [dummy,map] = nt2int('a'); %#ok
        maxval = max(map);
    end
if ischar(seq)
    try
        seq = nt2int(seq);
    catch
        result = false;
        return
    end
    
end

if acgtonly
    upperLimit = 4;
else
    upperLimit = maxval; % max(nt2int('ACGTUNRYKMSWBDHV*'));
end
result = ~any(any(seq <= 0 | seq > upperLimit | seq ~= floor(seq)));

end


