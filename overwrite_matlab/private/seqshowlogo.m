function seqshowlogo(varargin)
%SEQSHOWLOGO displays a Java seqlogo frame in a figure window

% Copyright 2004 The MathWorks, Inc.

isAA = false;
seqType = 'NT';
filename = 'seqlogo.png'; %#ok!
saveLogo = false;%#ok!
wtMatrix = [];
symbols = [];
% logoCell = {};

if nargin == 1 % Only pass in cell array and seqType default to NT
    logoCell = varargin{1};
    wtMatrix = logoCell{1,2};
    symbols = logoCell{1,1};
elseif nargin == 2 % Pass in cell array and isAA
    logoCell = varargin{1};
    isAA = varargin{2};
    wtMatrix = logoCell{1,2};
    symbols = logoCell{1,1};
elseif nargin == 3 % Pass in weight Matrix, list of symbols and isAA
    wtMatrix = varargin{1};
    symbols = varargin{2};
    isAA = varargin{3};
elseif nargin == 4 % Pass in weight Matrix, list of symbols, isAA and filename
    saveLogo = true;%#ok!
    wtMatrix = varargin{1};
    symbols = varargin{2};
    isAA = varargin{3};
    filename = varargin{4};%#ok!
end

if isAA
    seqType = 'AA';
end

% Get toolbox/matlab/icon path
iconPath = [toolboxdir('matlab'), filesep, 'icons', filesep];
import com.mathworks.toolbox.bioinfo.sequence.*;
import com.mathworks.mwswing.MJScrollPane;
import java.awt.Dimension;
% Create the viewer
b = SequenceViewer(wtMatrix, symbols, seqType);
awtinvoke(b,'addSeqLogo()');
scrollpanel = MJScrollPane(b, MJScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, MJScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

% Create a figure with the seqlogo panel on it and a uitoolbar
hFigure = figure( ...
            'WindowStyle', 'normal', ...
            'Menubar', 'none', ...
            'Toolbar', 'none', ...
            'Resize', 'on', ...
            'HandleVisibility', 'callback',...
            'IntegerHandle', 'off',...
            'NumberTitle','off',...
            'Tag', 'seqlogo',...
            'Name', 'Sequence Logo',...
            'DeleteFcn', @onLogoClosing);
        
% Set the figure widow size to fit the scrollPane
d = awtinvoke(scrollpanel, 'getPreferredSize()');
pos = getpixelposition(hFigure);
pos(3) = d.getWidth;
pos(4) = d.getHeight;
setpixelposition(hFigure,pos);
figurePosition = get(hFigure, 'Position');
[logoP, logoC] = javacomponent(scrollpanel, ...
    [0, 0, figurePosition(3), figurePosition(4)], ...
    hFigure);%#ok

set(logoC, 'units', 'normalized');

%  Add custom toolbar 
tb = uitoolbar(hFigure);
cicon= load([iconPath,'printdoc.mat']); % load cdata of print icon from toolbox/matlab/icon/printdoc.mat
a1=uipushtool(tb, ...
    'CData', cicon.cdata, ...
    'ClickedCallback', @printHandler); %#ok
cicon=load([iconPath,'savedoc.mat']);  % load cdata of save icon from the mat file
uipushtool(tb, ...
    'CData', cicon.cdata, ...
    'ClickedCallback', @saveHandler);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin nested functions 
function printHandler(varargin) 
% %      b.logoPrint;
awtinvoke(b, 'logoPrint()');
end

function saveHandler(varargin) 
     b.saveLogoDialog(wtMatrix, symbols, seqType);
% %      awtinvoke(b, 'saveLogoDialog([[D[CLjava/lang/String;)',  wtMatrix, symbols, seqType);
end

function onLogoClosing(hfig, event)%#ok
    if ~isempty(b)
        awtinvoke(b, 'cleanup()');
        delete(logoC);
        b=[];
    end
end

end %end of showlogo