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
%   $Revision: 1.10.6.10 $  $Date: 2010/12/22 16:19:23 $

acgtonly = false;
if nargin > 1
    if rem(nargin,2)== 0
        error(message('bioinfo:isnt:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'acgtuonly'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:isnt:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:isnt:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % others forces everything except ACTG to be the unknown value
                    if islogical(pval)
                        acgtonly = pval;
                    else
                        error(message('bioinfo:isnt:InvalidACGTOnly'));
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
