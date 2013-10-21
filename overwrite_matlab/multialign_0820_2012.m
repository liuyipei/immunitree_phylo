function iseqs = multialign(iseqs,tr,varargin)
%MULTIALIGN Progressive multiple sequence alignment.
%  (Yi's modified version)
%    I have added small noise to the position specific gap penalties
%
%  MA = MULTIALIGN(SEQS) performs progressive multiple alignment of the
%  sequences in SEQS. Pairwise distances between sequences are computed
%  after pairwise alignment with the Gonnet scoring matrix and then by
%  counting the proportion of sites at which each pair of sequences are
%  different (ignoring gaps). The guide tree is calculated by the
%  neighbor-joining method assuming equal variance and independence of
%  evolutionary distance estimates. SEQS is a vector of structures with the
%  fields 'Sequence' for the residues and 'Header' or 'Name' for the
%  labels. MA is the same structure as SEQS but with the field 'Sequence'
%  updated with the alignment. SEQS may also be a cell array of strings or
%  a char array, in such case MA is a char array with the output alignment
%  following the same order as the input.
%
%  MULTIALIGN(SEQS,TREE) uses the tree in TREE as a guide for the
%  progressive alignment. SEQS should have the same order as the leaves in
%  the TREE or use the field 'Header' (or 'Name') to identify the
%  sequences. TREE is a phylogenetic tree calculated with either SEQLINKAGE
%  or SEQNEIGHJOIN functions. 
%
%  MULTIALIGN(...,'WEIGHTS',METHOD) selects the sequence weighting method.
%  Weights are used to emphasize highly divergent sequences by scaling 
%  the score matrix and gap penalties. Closer sequences receive small
%  weights. Options are:  
%     'THG' (default)  - Thompson-Higgins-Gibson method using the
%                        phylogenetic tree branch distances weighted by
%                        their thickness. 
%     'equal'          - Assigns same weight to every sequence.
%
%  MULTIALIGN(...,'SCORINGMATRIX',SM) defines the [MxM] scoring matrix SM
%  to be used for the progressive alignment. SM can be a [MxMxN] array with
%  N user-defined scoring matrices. Match and mismatch scores are
%  interpolated from the series of scoring matrices considering the
%  distances between the two profiles (or sequences) being aligned; the
%  first matrix corresponds to the smallest distance and the last matrix to
%  the largest distance. Intermediate distances are calculated using linear
%  interpolation. The default is the BLOSUM80 to BLOSUM30 series for amino
%  acids or a fixed matrix NUC44 for nucleotides. When passing your own
%  series of scoring matrices make sure that all of them share the same
%  scale. SM may also be a cell array of strings with matrix names.
%
%  MULTIALIGN(...'SMINTERP',false) Turns off the linear interpolation of
%  the scoring matrices. Instead, each supplied scoring matrix is assigned
%  to a fixed range, depending on the distances between the two profiles
%  (or sequences) being aligned. Default is true.
%
%  MULTIALIGN(...,'GAPOPEN',P) defines the initial penalty for opening a
%  gap. P can be a scalar or a function specified using @. MULTIALIGN
%  passes four values to the function: the average score for two matched
%  residues (sm), the average score for two mismatched residues (sx), and,
%  the length of both profiles or sequences (len1 and len2). P defaults to
%  @(sm,sx,len1,len2) 5*sm. 
%
%  MULTIALIGN(...,'EXTENDGAP',Q) defines the initial penalty for extending
%  a gap. Q can be a scalar or a function specified using @. MULTIALIGN
%  passes four values to the function: the average score for two matched
%  residues (sm), the average score for two mismatched residues (sx), and
%  the length of both profiles or sequences (len1 and len2). Q defaults to
%  @(sm,sx,len1,len2) sm/4.  
%
%  MULTIALIGN(...,'DELAYCUTOFF',D) defines a threshold to delay the
%  alignment of divergent sequences whose closest neighbor is farther than
%  D * (median patristic distance between sequences). D defaults to the
%  unity, meaning that all sequences in which the closest sequence is
%  farther than the median distance are delayed.
%
%  MULTIALIGN(...,'USEPARALLEL',true) computes the pairwise alignments using
%  parfor loops.  If Parallel Computing Toolbox is installed and a
%  matlabpool is open, computation occurs in parallel. If Parallel
%  Computing Toolbox is not installed, or a matlabpool is not open,
%  computation uses parfor loops in serial mode. Defaults to false, or
%  serial computation with for loops. 
%
%  MULTIALIGN(...,'VERBOSE',true) turns on verbosity. Default is false.
%
%  The remaining input parameters are analogous to PROFALIGN and they are
%  used through every step of the progressive alignment of profiles. Get
%  help on PROFALIGN for more info. 
%
%   MULTIALIGN(...,'EXISTINGGAPADJUST',LOGICAL) Default is true.
%   MULTIALIGN(...,'TERMINALGAPADJUST',LOGICAL) Default is false.
%
%  Examples: 
%
%  % Align seven cellular tumor antigen p53 sequences
%  p53 = fastaread('p53samples.txt')
%  ma = multialign(p53,'verbose',true)
%  multialignviewer(ma)
% 
%  % Use an UPGMA phylogenetic tree instead as a guiding tree:
%  dist = seqpdist(p53,'ScoringMatrix',gonnet);
%  tree = seqlinkage(dist,'UPGMA',p53)
%  % and score the progressive alignment with the PAM family
%  ma = multialign(p53,tree,'ScoringMatrix',{'pam150','pam200','pam250'})
%  multialignviewer(ma)
%
%  % Promote terminations with gaps in the alignment
%  seqs = {'CACGTAACATCTC','ACGACGTAACATCTTCT','AAACGTAACATCTCGC'};
%  multialign(seqs,'terminalGapAdjust',true)
%             
%  See also BIRDFLUDEMO, HMMPROFALIGN, LOCALALIGN, MULTIALIGNREAD,
%  MULTIALIGNVIEWER, MULTIALIGNWRITE, NWALIGN, PROFALIGN, SEQCONSENSUS,
%  SEQNEIGHJOIN, SEQPROFILE, SHOWALIGNMENT. 

% References:
%   J.D. Thompson, D.G. Higgins, and T.J. Gibson. Nucleic Acids Res. (1994)
%   22(22):4673-4680.
%   R. Durbin, S. Eddy, A. Krogh, and G. Mitchison. Biological Sequence
%   Analysis. Cambridge UP, 1998.
%
% Undocumented constants to score the existing gap-gap and gap-residue
% cases that exist already in the profiles. These input arguments are
% analogous to the same undocumented input arguments in PROFALIGN. 
%
%  MULTIALIGN(...,'GAPGAPSCORE',GGS) default is @(sm,sx,len1,len2) 0.1*sm.
%  MULTIALIGN(...,'GAPRESSCORE',GRS) default is @(sm,sx,len1,len2) 0.1*sx.
%
%  OBSOLETE OPTIONS:
%
%  MULTIALIGN(...,'-v2',true) default is false. Use '-v2' to reproduce
%  results for earlier versions of Bioinformatics Toolbox 3.0.

%   Copyright 2005-2010 The MathWorks, Inc.
%   $Revision: 1.1.8.16 $  $Date: 2011/05/20 04:16:42 $


%-% defaults
useParallel = false;
dispLog = false;
profalignArgs = {};
ogIsFunctionHandle = true;
egIsFunctionHandle = true;
ogF = @(sm,sx,l1,l2) 5*sm;
egF = @(sm,sx,l1,l2) sm/4;
distSeqsToDelay = 1;
weightMethod = 'thg'; % Thompson-Higgins-Gibson
interpolateSM = true;
seqNamesGiven = true;
defaultSM = true;
calculateTree = false;

%-% check input sequences 
% they can be an array of chars, a vector of string cells, or, a vector
% array of structures with the field 'Sequence'
if ischar(iseqs)
    seqs = mat2cell(iseqs,ones(1,size(iseqs,1)),size(iseqs,2));
elseif isstruct(iseqs) && isfield(iseqs,'Sequence')
    seqs = {iseqs(:).Sequence};
elseif isstruct(iseqs) && isfield(iseqs,'sequence')  % other functions also allow sequence
    seqs = {iseqs(:).sequence};
    iseqs = rmfield(iseqs,'sequence');
    warning('bioinfo:multialign:StructSeqCapitalization',...
            ['Field names in MATLAB structures are case sensitive.\nThe input ',...
            'structure contains a ''sequence'' field. However, most \nfunctions in the ',...
            'toolbox create structures with a ''Sequence'' field.\nUsing mixed ',...
            'capitalization can lead to unexpected results.']);
elseif iscellstr(iseqs)
    seqs = iseqs;
else
    error(message('bioinfo:multialign:IncorrectInputType'))
end
seqs = strtrim(seqs(:));     % trim sequences
seqs = strrep(seqs,' ','-'); % align chars can only be '-'
seqs = strrep(seqs,'.','-'); % align chars can only be '-'
numSeqs = numel(seqs);
if numSeqs<3 
    error('bioinfo:multialign:InvalidNumberOfSequences',...
          'At least 3 input sequences must be supplied to MULTIALIGN.')
end

% if the input looks like a multiple alignment, first remove all columns
% which only have gaps
if numel(unique(cellfun('length',seqs)))==1
    gapCols = find(all(char(seqs)=='-'));
    if numel(gapCols)>0
        for i = 1:numSeqs
            seqs{i}(gapCols) = []; 
        end
    end
end

% try to guess which type of sequences we have
if isnt(strcat(seqs{:}))
    isAminoAcid = false;
    alpha = 'nt';
elseif isaa(strcat(seqs{:}))
    isAminoAcid = true;
    alpha = 'aa';
else
    error(message('bioinfo:multialign:InvalidInputOfSequences'))
end

%-% get the original names of each sequence
if isfield(iseqs,'Header') 
    seqNames = {iseqs(:).Header};
elseif isfield(iseqs,'Name')  
    seqNames = {iseqs(:).Name};
elseif isfield(iseqs,'LocusName')  
    seqNames = {iseqs(:).LocusName};
else 
    seqNames = strread(sprintf('seq %d\n',1:numSeqs),'%s','delimiter','');
    seqNamesGiven = false;
end

%-% check if the second input is a phylogenetic tree
if nargin==1 || ~strcmp(class(tr),'phytree')
    calculateTree = true;
    if nargin>1
        varargin = {tr,varargin{:}};
    end
end
   
%-% Check optional PV input arguments
nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2) == 1
        error(message('bioinfo:multialign:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'scoringmatrix','gapopen','extendgap','verbose',...
              'existinggapadjust','terminalgapadjust','delaycutoff',...
              'weights','sminterp','gapgapscore','gapresscore','useparallel',...
              '-v2'}; 
              
    for j=1:2:nvarargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:multialign:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:multialign:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % scoring matrix
                    if isnumeric(pval)
                        SM = pval;
                        defaultSM = false;
                    else
                        if ischar(pval)
                            pval = {pval};
                        end
                        if iscellstr(pval)
                            n = numel(pval);
                            try
                                for i = 1:n
                                    [SMt,SMi] = feval(lower(pval{i}));
                                    if isAminoAcid
                                        SM(:,:,i) = SMt(1:20,1:20) * SMi.Scale;
                                    else
                                        SM(:,:,i) = SMt(1:4,1:4) * SMi.Scale;
                                    end
                                end
                            catch allExceptions
                                error(message('bioinfo:multialign:InvalidStringScoringMatrix', pval{ i }));
                            end
                        else
                            error(message('bioinfo:multialign:InvalidScoringMatrix'));
                        end
                        defaultSM = false;
                    end
                case 2 %gap open penalty
                    if isscalar(pval) && isnumeric(pval)
                        ogIsFunctionHandle = false;
                        ogS = pval;
                    elseif isa (pval, 'function_handle')
                        ogF = pval;
                    else
                        error(message('bioinfo:multialign:InvalidGapOpen'))
                    end
                case 3 %extend gap penalty
                    if isscalar(pval) && isnumeric(pval)
                        egIsFunctionHandle = false;
                        egS = pval;
                    elseif isa (pval, 'function_handle')
                        egF = pval;
                    else
                        error(message('bioinfo:multialign:InvalidGapExtend'))
                    end
                case 4 % verbose
                    dispLog = opttf(pval);
                    if isempty(dispLog)
                        error(message('bioinfo:multialign:showscoreInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 5 % existinggapadjust
                    if isempty(opttf(pval))
                        error(message('bioinfo:multialign:existingGapAdjustInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                    profalignArgs(end+1:end+2)={'existingGapAdjust',pval};
                case 6 % terminalgapadjust  
                    if isempty(opttf(pval))
                        error(message('bioinfo:multialign:endGapAdjustInputOptionNotLogical', upper( char( okargs( k ) ) )));
                    end
                    profalignArgs(end+1:end+2)={'terminalGapAdjust',pval};
                case 7 % delay cutoff
                    if isscalar(pval) && pval>0
                        distSeqsToDelay = pval;
                    else
                        error(message('bioinfo:multialign:invalidDelayCutoff'))
                    end
                case 8 % weights
                    weightMethods = {'thg','equal'};
                    weightMethod = find(strncmpi(pval,weightMethods,numel(pval)));
                    if isempty(weightMethod) 
                        error(message('bioinfo:multialign:NotValidWeighMethod'))
                    end
                    weightMethod = weightMethods{weightMethod};
                case 9 % interpolte scoring matices 
                    interpolateSM = opttf(pval);
                    if isempty(interpolateSM)
                        error(message('bioinfo:multialign:interpolateSMOptionNotLogical', upper( char( okargs( k ) ) )));
                    end  
                case 10 % gapgapscore
                    if isscalar(pval) && isnumeric(pval)
                        profalignArgs(end+1:end+2)={'gapgapscore',pval};
                    elseif isa (pval, 'function_handle')
                        profalignArgs(end+1:end+2)={'gapgapscore',pval};
                    else
                        error(message('bioinfo:multialign:InvalidGapGapScore'))
                    end
                case 11 % gapresscore
                    if isscalar(pval) && isnumeric(pval)
                        profalignArgs(end+1:end+2)={'gapresscore',pval};
                    elseif isa (pval, 'function_handle')
                        profalignArgs(end+1:end+2)={'gapresscore',pval};
                    else
                        error(message('bioinfo:multialign:InvalidGapResScore'))
                    end
                case 12 % useParallel
                    useParallel = opttf(pval);
                    if isempty(useParallel)
                        error(message('bioinfo:multialign:useParallelNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 13 % -v2
                    error('bioinfo:multialign:Obsolete',...
                    'Input argument ''-v2'' is obsolete.');
%                    profalignArgs(end+1:end+2)={'-v2',pval};   
            end
        end
    end
end


%------------------------- get the guiding tree ---------------------------
% if tree was not given calculate it, otherwise just check consistency
if calculateTree 
    if numel(unique(seqNames))~=numSeqs
        seqNames = strread(sprintf('seq %d\n',1:numSeqs),'%s','delimiter','');
        warning(message('bioinfo:multialign:NotUniqueSequenceNames'))
    end

    d = seqpdist(seqs,'indels','pairwise-delete','method','p-dist',...
        'gapopen',10,'extend',1,'ScoringMatrix',gonnet,...
        'alphabet',alpha,'UseParallel',useParallel);
                  
    tr = seqneighjoin(d,'equivar',seqNames);
    [trPtrs trNames] = get(tr,'pointers','leafnames');
    hOrder = seqmatch(trNames,seqNames,'exact',true);
    % re-order tree pointer to have the same indices as the input sequences
    trPtrs(trPtrs<=numSeqs) = hOrder(trPtrs(trPtrs<=numSeqs));
else
    [trPtrs numLeaves trNames] = get(tr,'pointers','numleaves','leafnames');
    if numSeqs~=numLeaves
        error('bioinfo:multialign:InvalidNumberOfSequences',...
              'The number of sequences differs from the number of leaves in the TREE.')
    end
    if seqNamesGiven
        warnSt = warning('off','bioinfo:seqmatch:StringNotFound');
        hOrder = seqmatch(trNames,seqNames,'exact',true);
        warning(warnSt)
        % check consistency of the tree leaf names and the sequence names
        if ~isequal(1:numSeqs,sort(hOrder)')
            error(message('bioinfo:multialign:InvalidNamesOfSequences'))
        end
        % re-order tree pointer to have the same indices as the input
        % sequences
        trPtrs(trPtrs<=numSeqs) = hOrder(trPtrs(trPtrs<=numSeqs));
    else
        % in this case the order of the input sequences should have been the
        % same as the order of the tree leaves, since seqNames were not
        % given we cannot verify consistency
        hOrder = (1:numSeqs)';
    end
end

%----------------------- Prepare some variables ---------------------------

% prepare the default scoring matrices (if not given at the input)
if defaultSM
    if isAminoAcid
        SM = zeros(20,20,11);
        j=1;
        for i=80:-5:30
            [t ti] = blosum(i);
            SM(:,:,j) = t(1:20,1:20)*ti.Scale;
            j=j+1;
        end
    else
        [SM,SMi] = nuc44;
        SM = SM * SMi.Scale;
    end
end


% get tree distances
trDist = pdist(tr,'squareform',true);
trDist(hOrder,hOrder) = trDist; % reorder as input sequences 
trDistMax = max(trDist(:));
trDistAvg = median(trDist(~eye(numSeqs)));
trDist = trDist + diag(inf(numSeqs,1));
trDistMin = min(trDist(:));

% find distance to closest neighbor for each sequence
toDelay = min(trDist);
toDelay = toDelay>max(min(toDelay),trDistAvg*distSeqsToDelay);
     % the max is used to ensure that not all sequences are delayed, at
     % least the two closest will be aligned by the tree
toDelay(numSeqs*2-1) = false; % do not delay branches

% seqPtrs are pointers to aid the construction of the final alignment
seqPtrs = cell(1,numSeqs);
for i = 1:numSeqs
    seqPtrs{i} = (1:numel(seqs{i}))';
end

% convert initial sequences to profiles and allocate space for the new
% profiles, grps indicates which sequences belong to each profile
profs = cell(1,numSeqs*2-1);
grps =  cell(1,numSeqs*2-1);
for i = 1:numSeqs
    profs{i} = seqprofile(seqs{i},'gaps','all','alphabet',alpha,'ambiguous','count');
    grps{i} = i;
end

% Weight sequences:
% now we only use the Thompson-Higgins-Gibson method.
if strcmp(weightMethod,'thg') && trDistMax>0
    W(hOrder) = weights(tr); % weights reordered as the input sequences
    for i = 1:numSeqs
        profs{i} = profs{i} .* W(i);
    end
end

% prepare and check scoring matrix
siz = size(SM);
if isAminoAcid
    if any(siz([1 2])<20) || siz(1)~=siz(2)
        error(message('bioinfo:multialign:invalidAAScoringMatrix'))
    else
        SM = SM(1:20,1:20,:);
    end
else
    if any(siz([1 2])<4) || any(siz([1 2])>15) || siz(1)~=siz(2)
        error(message('bioinfo:multialign:invalidNTScoringMatrix'))
    else
        SM = SM(1:4,1:4,:);
    end
end

% prepare vector to interpolate scoring matrices
nSM = size(SM,3);
if nSM > 1
    if trDistMax>trDistMin
        x = (0:nSM-1) * (trDistMax-trDistMin)/(nSM-1) + trDistMin;
    else % if all distances are equal then we assume that we have a highly
         % conserved alignment and create an x vector such that at the time of
         % interpolating it picks the first scoring matrix
         x = [trDistMin realmax./(nSM-1:-1:1)];
    end
end

%----------------------- Algorithm starts here ----------------------------
if nSM == 1  % there is only one scoring matrix, do not need to interp
    iSM = SM;
elseif isAminoAcid
    iSM = zeros(20,20); % temporal array with the correct size
    SM = reshape(SM,400,nSM)';
else
    iSM = zeros(4,4);   % temporal array with the correct size
    SM = reshape(SM,16,nSM)';
end

% do progressive alignment
for i = 1:numSeqs-1
    ptrs = trPtrs(i,:);
    if all(toDelay(ptrs))
        % join grps of sequences
        grps{i+numSeqs} = cell2mat(grps(ptrs));
        toDelay(i+numSeqs) = true;
        if dispLog  % echo partial alignments
            disp(sprintf('Branch/Sequences delayed: %s',...
                sprintf('%s ',seqNames{grps{i+numSeqs}})))
            disp(' ')
        end
    elseif any(toDelay(ptrs))
        % delay one branch, pass the other to the upper level
        profs(i+numSeqs) = profs(ptrs(~toDelay(ptrs)));
        grps(i+numSeqs) = grps(ptrs(~toDelay(ptrs)));
        if dispLog  % echo partial alignments
            disp(sprintf('Branch/Sequences delayed: %s',...
                sprintf('%s ',seqNames{grps{ptrs(toDelay(ptrs))}})))
            disp(' ')
        end
    else % align profiles of both branches
        % find an adjusted scoring matrix for the profile alignment
        d = trDist(grps{ptrs(1)},grps{ptrs(2)});
        d = median(d(:));
        if nSM>1
            if interpolateSM
                iSM(:) = interp1(x,SM,d,'linear');
            else
                iSM(:) = interp1(x,SM,d,'nearest');            
            end
        end
        % calculate some statistics that might change some parameters in
        % the profile alignment
        if ogIsFunctionHandle || egIsFunctionHandle 
            sm = mean(diag(iSM));
            sx = sum(sum(iSM-diag(diag(iSM))))/prod(size(iSM)-[0 1]);
            l1 = size(profs{ptrs(1)},2);
            l2 = size(profs{ptrs(2)},2);
        end
        % set initial gap values
        if ogIsFunctionHandle
            ogS = ogF(sm,sx,11,12);
        end
        if egIsFunctionHandle
            egS = egF(sm,sx,11,12);
        end
        
        %% Yi -- adding some regularization
        % ogS and egS are just scalars in multialign -- make them into
        % position specific vectors, acoording to spec: a 2-item cell
        % array, where each item is a vector of length greater than each of
        % the two profiles by one.
        ogS_cellArr = {ogS + (0:size(profs{ptrs(1)},2))*100*eps, ogS + (0:size(profs{ptrs(2)},2))*100*eps};
        egS_cellArr = {egS + (0:size(profs{ptrs(1)},2))*100*eps, egS + (0:size(profs{ptrs(2)},2))*100*eps};
        
        % profile alignment
        [profs{i+numSeqs} h1 h2] = profalign(profs{ptrs},...
                                  'gapOpen',ogS_cellArr,'extendGap',egS_cellArr,...
                                  'scoringMatrix',iSM,...
                                  profalignArgs{:});
        % join grps of sequences
        grps{i+numSeqs} = cell2mat(grps(ptrs));
        % adjust pointers of the sequences of the first profile
        for k = grps{ptrs(1)}
            seqPtrs{k}=h1(seqPtrs{k});
        end
        % adjust pointers of the sequences of the second profile
        for k = grps{ptrs(2)}
            seqPtrs{k}=h2(seqPtrs{k});
        end

        if dispLog  % echo partial alignments
            maTemp = repmat('-',numel(grps{i+numSeqs}),size(profs{i+numSeqs},2));
            for j= 1:numel(grps{i+numSeqs})
                maTemp(j,seqPtrs{grps{i+numSeqs}(j)})=seqs{grps{i+numSeqs}(j)};
            end
            disp(sprintf('Branch No %d',i))
            if ogIsFunctionHandle || egIsFunctionHandle
                disp(sprintf('Match scr avg: %0.4f  Mismatch scr avg: %0.4f  Profile lengths: %d %d',sm,sx,l1,l2))    
            end
                        
            disp(sprintf('Open gap: %0.4f   Gap extend: %0.4f   Tree distance: %0.4f',ogS,egS,d))                        
            disp(sprintf('Aligned sequences:'))
            disp([char(seqNames{grps{i+numSeqs}}) repmat('  ',numel(grps{i+numSeqs}),1) maTemp])
            disp(' ')
        end
    end
end

% align delayed sequences
rootInd = numSeqs*2-1;
for i = find(toDelay(1:numSeqs))
    % calculate the interpolated scoring matrix
    d = trDist(grps{ptrs(1)},grps{ptrs(2)});
    d = median(d(:));
    if nSM>1
        if interpolateSM
            iSM(:) = interp1(x,SM,d,'linear');  
        else
            iSM(:) = interp1(x,SM,d,'nearest');
        end
    end
    % set initial gap values
    if ogIsFunctionHandle || egIsFunctionHandle
        sm = mean(diag(iSM));
        sx = sum(sum(iSM-diag(diag(iSM))))/prod(size(iSM)-[0 1]);
        l1 = size(profs{ptrs(1)},2);
        l2 = size(profs{ptrs(2)},2);
    end
    if ogIsFunctionHandle
        ogS = ogF(sm,sx,11,12);
    end
    if egIsFunctionHandle
        egS = egF(sm,sx,11,12);
    end
    
    %% Yi -- adding some regularization
    % ogS and egS are just scalars in multialign -- make them into
    % position specific vectors, acoording to spec: a 2-item cell
    % array, where each item is a vector of length greater than each of
    % the two profiles by one.
    ogS_cellArr = {ogS + (0:size(profs{i},2))*100*eps, ogS + (0:size(profs{rootInd},2))*100*eps};
    egS_cellArr = {egS + (0:size(profs{i},2))*100*eps, egS + (0:size(profs{rootInd},2))*100*eps};
    
    [profs{rootInd} h1 h2] = profalign(profs{[i,rootInd]},...
        'gapOpen',ogS_cellArr,'extendGap',egS_cellArr,...
        'scoringMatrix',iSM,profalignArgs{:});
    
    for k = grps{rootInd}
        seqPtrs{k} = h2(seqPtrs{k});
    end
    seqPtrs{i} = h1(seqPtrs{i});
    grps{rootInd} = [grps{rootInd} i];
    if dispLog  % echo partial alignments
        maTemp = repmat('-',numel(grps{rootInd}),size(profs{rootInd},2));
        for j= 1:numel(grps{rootInd})
            maTemp(j,seqPtrs{grps{rootInd}(j)})=seqs{grps{rootInd}(j)};
        end
        disp(sprintf('Aligning delayed sequence %s',seqNames{grps{i}}))
        if ogIsFunctionHandle || egIsFunctionHandle
            disp(sprintf('Match scr avg: %0.4f  Mismatch scr avg: %0.4f  Profile lengths: %d %d',sm,sx,l1,l2))
        end
        disp(sprintf('Open gap: %0.4f   Gap extend: %0.4f   Tree distance: %0.4f',ogS,egS,d))
        disp([char(seqNames{grps{rootInd}}) repmat('  ',numel(grps{rootInd}),1) maTemp])
        disp(' ')
    end
end

%-% Build the output alignment
alignLength = size(profs{rootInd},2);
if isstruct(iseqs) % output is a structure
    for i = 1:numSeqs
        iseqs(i).Sequence = repmat('-',1,alignLength);
        iseqs(i).Sequence(seqPtrs{i}) = seqs{i};
    end
    gapCols = find(all(char(iseqs.Sequence)=='-'));
    if numel(gapCols)>0
        for i = 1:numSeqs
            iseqs(i).Sequence(gapCols) = [];
        end
    end
else               % output is a char array
    iseqs = repmat('-',numSeqs,alignLength);
    for i = 1:numSeqs
        iseqs(i,seqPtrs{i}) = seqs{i};
    end
    iseqs(:,all(iseqs=='-')) = [];
end
