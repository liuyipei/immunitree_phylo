function color = setcolorpref(pval, varargin)
% SETCOLORPREF Generic color options handler.
%
% C = SETCOLORPREF(PVAL,  DEFCOLOR) parse the color input and returns color
% C in hex format to be used in html browser or in Java applications. If
% PVAL is not a right format. Throw a warning and return a default color
% DEFCOLOR or black if not specified.

%   Copyright 2003-2009 The MathWorks, Inc.
%   $Revision: 1.4.4.5 $  $Date: 2010/12/22 16:19:30 $

if nargin > 1
    defColor = varargin{1};
else
    defColor = '000000';
end
color = '';
if isnumeric(pval) && isvector(pval) && max(size(pval))== 3
    if max(pval) <= 1 && min(pval)>=0
        pval = fix(pval*255);
    end
    color = sprintf('%02x%02x%02x',pval(1), pval(2),pval(3));
elseif ischar(pval)
    switch lower(pval)
        case 'r'
            color = 'FF0000';
        case 'g'
            color = '00CC00';
        case 'b'
            color = '0000FF';
        case 'm'
            color = 'FF00FF';
        case 'c'
            color = '00FFFF';
        case 'y'
            color = 'CCCC00';
        case 'w'
            color = 'FFFFFF';
        case 'k'
            color = '000000';
        case 'red'
            color = 'FF0000';
        case 'green'
            color = '00CC00';
        case 'blue'
            color = '0000FF';
        case 'magenta'
            color = 'FF00FF';
        case 'cyan'
            color = '00FFFF';
        case 'yellow'
            color = 'CCCC00';
        case 'white'
            color = 'FFFFFF';
        case 'black'
            color = '000000';
    end
elseif numel(pval) == 6  % or HTML style string
        pval = upper(pval);
        if isempty(regexp(pval,'([^A-F0-9])','once'))
            color = pval;
        end
% % else
% %     bioinfoprivate.bioerror(mfilename,...
% %         'InvalidColorInput',...
% %         'COLOR must be a three-element RGB vector or one of the predefined names.');
end

if isempty(color)
    color = defColor;
    warning(message('bioinfo:setcolorpref:BadColorSpec'));
end
