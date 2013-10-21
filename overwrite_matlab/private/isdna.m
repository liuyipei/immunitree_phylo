function result = isdna(seq,varargin)
%ISDNA True for DNA sequences.
%   ISDNA(SEQ) returns 1 for a DNA sequence, 0 otherwise. Valid symbols are
%   A,C,G,T,N,R,Y,K,M,S,W,B,D,H,V and *.
%
%   ISDNA(...,'ACGTOnly',true) returns 1 only if the sequence contains
%   A,C,G and T only.   
%
%   See also ISNT, ISRNA, ISAA.

%   Copyright 2002-2004 The MathWorks, Inc.
%   $Revision: 1.5.4.2 $  $Date: 2005/06/09 21:57:24 $

if ischar(seq)
    if any(lower(seq) == 'u')
        result = false;
        return
    end
end
result = isnt(seq,varargin{:});
