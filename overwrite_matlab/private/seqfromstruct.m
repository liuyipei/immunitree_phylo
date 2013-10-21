function seq = seqfromstruct(seqStruct)
%SEQFROMSTRUCT returns Sequence field of a structure.
%
%   SEQFROMSTRUCT(STRUCTURE) returns the field Sequence from STRUCTURE. If
%   no Sequence field exists but a field with some other capitalization of
%   Sequence does exist, then this field is returned but a warning is
%   given.
%
%   If the input contains a nested structure then the leaf node is
%   returned.
%
%   The function gives an MException if it encounters an array of structures.

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2010/12/22 16:19:28 $

if numel(seqStruct) > 1
    xcptn = MException('bioinfo:seqfromstruct:SeqStructArray',...
        ['The input structure contains multiple sequences.\n',...
        'Please specify a single sequence.']);
    xcptn.throwAsCaller;
end
% if the struct has a field Sequence then we use it
try
    seq = seqStruct.Sequence;
catch
    % if not we check for fields with different capitalization
    fields = fieldnames(seqStruct);
    matches = find(strcmpi(fields,'sequence'));
    if numel(matches) == 1
        seq = seqStruct.(fields{matches});
        warning('bioinfo:seqfromstruct:StructSeqCapitalization',...
            ['Field names in MATLAB structures are case sensitive.\nThe input ',...
            'structure contains a ''%s'' field. However, most \nfunctions in the ',...
            'toolbox create structures with a ''Sequence'' field.\nUsing mixed ',...
            'capitalization can lead to unexpected results.'], fields{matches});
    elseif numel(matches) > 1
        xcptn = MException('bioinfo:seqfromstruct:StructMultiSeqFields',...
            ['The input structure contained multiple sequence fields.\n',...
            'Please specify a single sequence.']);
        xcptn.throwAsCaller;
    else
        xcptn = MException('bioinfo:seqfromstruct:StructNoSequence',...
            'Input value is a structure that does not contain a Sequence field.');
        xcptn.throwAsCaller;
    end
end
if isstruct(seq)
    % call recursively to find the sequence
    try
        seq = seqfromstruct(seq);
    catch Mxcptn
        if strcmp(Mxcptn.identifier,'bioinfo:seqfromstruct:StructNoSequence')
            xcptn = MException('bioinfo:seqfromstruct:BadSequenceField',...
                'The Sequence field in the input structure contains a sub-structure without a Sequence field.');
            xcptn.throwAsCaller;
        end
        Mxcptn.throwAsCaller;
    end
end