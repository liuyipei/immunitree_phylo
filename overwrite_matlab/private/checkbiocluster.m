function checkbiocluster(varargin)
% CHECKBIOCLUSTER check the existence and availability of a cluster.
%
%   THIS FUNCTION IS OBSOLETE. USE A MATLABPOOL INSTEAD FOR PARALLEL OR
%   DISTRIBUTED PROCESSING.

%   Copyright 2003-2010 The MathWorks, Inc.
%   $Revision: 1.1.8.4 $  $Date: 2010/12/22 16:19:16 $

error(message('bioinfo:checkbiocluster:Obsolete'));

% if ~isa(jm,'distcomp.jobmanager')
%     if isempty(ver('distcomp'))
%         error('bioinfo:checkbiocluster:NoDCT',...
%               'The Parallel Computing Toolbox is required to distribute the algorithm into a cluster of computers.')
%     else
%         error('bioinfo:checkbiocluster:InvalidJobManager',...
%             'Invalid job manager. Use findResource from the Parallel Computing Toolbox to find a valid job manager.')
%     end
% end
% state = jm.State;
% if strcmp(state,'unavailable')
%     error('bioinfo:checkbiocluster:JobManagerUnavailable',...
%         'Job manager is unavailable.')
% end
% if ~strcmp(state,'running')
%     error('bioinfo:checkbiocluster:JobManagerNotRunning',...
%         'Job manager ''%s'' state is not ''running''.',jm.Name)
% end
% 
% u = get(jm.Jobs,'State');
% 
% if ~waitInQueue && (any(strcmp(u,'queued')) || any(strcmp(u,'running')))
%     error('bioinfo:checkbiocluster:NoIddleWorkers',...
%           'Job manager ''%s'' does not have available resources. To wait in the queue set the option WAITINQUEUE to true.',jm.Name)
% end
% 

