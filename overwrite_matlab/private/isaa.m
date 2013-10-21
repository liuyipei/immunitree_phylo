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

