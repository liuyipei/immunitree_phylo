function hg_handles = dist_node_draw(node_handle)
%CUSTOMNODEDRAW User function example to draw customized nodes.
%
%  Utilize a user specified function to customize the node drawing of a 
%  BIOGRAPH object. The function must have the following form:
%
%     function  HG_HANDLES = CUSTOMNODEDRAW(NODE_HANDLE)
%
%  taking as input argument a @biograph/@node object. The function
%  CUSTOMNODEDRAW uses the @biograph/@node properties for drawing the node,
%  such as the "Position", the "Size" and the "Label" of the node. The user
%  can store application data into the property "UserData" before creating
%  the graph layout, then CUSTOMNODEDRAW can utilize it to enhance the node
%  rendered information, for example confidence values, rates or
%  probability distributions.
%
%  The handles of all objects drawn by CUSTOMNODEDRAW must be returned into
%  the vector HG_HANDLES, so the @biograph object deletes them when
%  necessary.


%   Copyright 2003-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2008/03/28 15:21:03 $


% Get the axes handle where all hg objects are attached, all hg objects
% should specifically set "haxes" as the parent to ensure that the
% customized node is rendered in the proper figure.
haxes = node_handle.up.hgAxes;

% Get the current scale between the computed layout and rendered graph
biographScale = node_handle.up.Scale;

% Get the size of the node
rx = node_handle.Size(1);
ry = node_handle.Size(2);

% Get the center of the node
x = node_handle.Position(1)*biographScale;
y = node_handle.Position(2)*biographScale;

% In this example you create a pie chart with the vector stored in
% node_handle.UserData.Distribution

if isempty(node_handle.UserData) || ...
        ~isfield(node_handle.UserData,'Distribution') || ...
        isempty(node_handle.UserData.Distribution)

    % could not find data to draw the pie, then just draw a circle
    t = [0:pi/25:2*pi];
    px = [rx/2*cos(t)+x];
    py = [ry/2*sin(t)+y];
    hg_handles = zeros(2,1); % initialize output handles
    hg_handles(1) = patch(px,py,[1 0 0],'Parent',haxes);
    hg_num = 1;
    
else % draw the pie chart
     % 2012/07: Yi: introduce matrix format 
     % (compatible with Joni's original one-row format)
     % In new format, sums of columns correspond to large outer-slices
     % each row subdivies the outer slice into inner slices
    
    d_coarse = sum(node_handle.UserData.Distribution, 1);    
    % turn into a column, if necessary, and normalize    
    d_coarse = d_coarse(:)./sum(d_coarse(:)); 

    alpha_coarse = [0;cumsum(d_coarse)]*2*pi;

    % Draw each slice of the pie chart
    num_slices = numel(d_coarse);
    num_subslices_per_slice = size(node_handle.UserData.Distribution, 1);

    colors = jet(num_slices);    
    % Joni/Yi 01/05/2012
    % added by yoni so that the color of the first element will be red
    %if num_slices <= 4
    %    colors = fliplr(jet(4));
    %else
    if num_slices <= 8
        colors = jet(8);
    end
    
    hg_handles = zeros(2*num_slices+1,1); % initialize output handles
    hg_num = 0;

    for i = 1:num_slices
        if d_coarse(i) == 0 || ~isfinite(d_coarse(i))            
            px = [x x x];
            py = [y y y];
            continue;
        else
            t = [alpha_coarse(i):pi/12:alpha_coarse(i+1) alpha_coarse(i+1)];
            if num_subslices_per_slice > 1 % leave room for inner pie
                trev = t(end:-1:1);            
                px = [rx/2*cos(t)+x rx/3*cos(trev)+x];
                py = [ry/2*sin(t)+y ry/3*sin(trev)+y];
                px = [px px(1)];
                py = [py py(1)];
            else % just do the regular pie slice
                px = [x rx/2*cos(t)+x x];
                py = [y ry/2*sin(t)+y y];    
            end
            hg_handles(hg_num+1) = patch(px,py,colors(i,:),'Parent',haxes);
            hg_num = hg_num+1;
        end                
    end
    
    % create subslices
    if num_subslices_per_slice > 1
        num_subslices = numel(node_handle.UserData.Distribution);
        % go down first col, then second, then etc
        d_fine = node_handle.UserData.Distribution(:)./sum(node_handle.UserData.Distribution(:)); 
        subslice_colors = repmat(jet(num_subslices_per_slice),num_slices,1);
        alpha_fine = [0;cumsum(d_fine)]*2*pi;
        for i = 1:num_subslices
            if d_fine(i) == 0 || ~isfinite(d_fine(i))            
                px = [x x x];
                py = [y y y];
                continue;
            else
                t = [alpha_fine(i):pi/25:alpha_fine(i+1) alpha_fine(i+1)];
                px = [x rx/3*cos(t)+x x];
                py = [y ry/3*sin(t)+y y];                
                hg_handles(hg_num+1) = patch(px,py,subslice_colors(i,:),'Parent',haxes);        
                hg_num = hg_num+1;
            end
        end
    end

end

% Place the label in the lower-right corner of the node
hg_handles(hg_num+1) = text(x+rx/2,y-ry/2,node_handle.Label,...
    'HorizontalAlignment','Left',...
    'VerticalAlignment','Middle',...
    'Fontsize', node_handle.FontSize, ...
    'Interpreter','none',...
    'Parent',haxes);
hg_num = hg_num+1;
hg_handles = hg_handles(1:hg_num);