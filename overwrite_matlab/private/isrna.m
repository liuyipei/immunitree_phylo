function result = isrna(seq,varargin)
%ISRNA True for RNA sequences.
%   ISRNA(SEQ) returns 1 for a RNA sequence, 0 otherwise. Valid symbols are
%   A,C,G,U,N,R,Y,K,M,S,W,B,D,H,V and *.
%
%   ISRNA(...,'ACGTOnly',true) returns 1 only if the sequence contains
%   A,C,G and U only.   
%
%   See also ISNT, ISDNA, ISAA.

%   Copyright 2002-2004 The MathWorks, Inc.
%   $Revision: 1.5.4.2 $  $Date: 2005/06/09 21:57:26 $

if ischar(seq)
    if any(lower(seq) == 't')
        result = false;
        return
    end
end
result = isnt(seq,varargin{:});
