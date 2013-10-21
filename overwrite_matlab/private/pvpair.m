function [k, pval] = pvpair(pname, theVal, okargs,mfile)
% PVPAIR Helper function that looks for partial matches of parameter names
% in a list of inputs and returns the parameter/value pair and matching
% number.
%
% [K, PVAL] = PVPAIR(PNAME, THEVAL, OKARGS) given input string PNAME,
% and corresponding value, THEVAL, finds matching name in the OKARGS list.
% Returns K, the index of the match, and PVAL, the parameter value.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.2 $   $Date: 2010/12/22 16:19:26 $

k = find(strncmpi(pname, okargs,numel(pname)));
if numel(k) == 1
    pval = theVal;
    return
end

if isempty(k)
    xcptn = MException(sprintf('bioinfo:%s:UnknownParameterName',mfile),...
        'Unknown parameter name: %s.',pname);
    xcptn.throwAsCaller;

elseif length(k)>1
    xcptn = MException(sprintf('bioinfo:%s:AmbiguousParameterName',mfile),...
        'Ambiguous parameter name: %s.',pname);
    xcptn.throwAsCaller;
end
