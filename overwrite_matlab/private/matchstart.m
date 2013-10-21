function tf = matchstart(string,pattern)
%MATCHSTART matches start of string with pattern, ignoring spaces

% Copyright 2003-2004 The MathWorks, Inc.
% $Revision: 1.1.12.1 $   $Date: 2004/12/24 20:42:38 $

tf = ~isempty(regexp(string,['^(\s)*?',pattern],'once'));
