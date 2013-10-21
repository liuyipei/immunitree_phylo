% Phylo version

function h = view_tree(T, sequences, ...    
    pretty_max_colors, time_stamp_data, Z,...
    germline, sublabels)
    h = [];
    if ~exist('germline','var')
        germline.frame_shift = 0;
        germline.seq = 5*ones(1, size(sequences,2));
	    germline.mutcount_end_indx = size(sequences,2);
    end
    if ~exist('sublabels', 'var')
        sublabels = [];
    end
    
    % mark deleted nodes in the sequences array as well
    ix = find(T(:,1) == -1);
    sequences(ix,:) = -1;
    
%     % ascii print the tree
    [X nodes depth] = order_tree(T);
    labels = cell(size(X,1), 1);
    ids = cell(size(X,1), 1);

    % Store the mutation strings in the Description of the Biograph nodes
    [~, mut, ~] = ...
      annotate_mutations_on_tree(T, sequences, true, ...
        germline.mutcount_end_indx, false);
    do_pairwise_mutstr = false
    [nMuts, ~, mut_str] = ...
      annotate_mutations_on_tree(T, sequences, true, ...
        germline.mutcount_end_indx, do_pairwise_mutstr);

    mut_str

    for i=1:size(X,1)
        if X(i,1) == 0, continue; end
        k = depth(i);
        t = nodes(i);
        mut_ct_str = sprintf('%d mutations', nMuts(t));
        labels{t} = sprintf('%d/%d', T(t,2), T(t,3));
        ids{t} = sprintf('%d', t);
        for j=1:k-1, fprintf(' \t'); end
        fprintf('%s (%s)  %s\n', ids{t}, labels{t}, mut_ct_str);
    end
    nMuts(nMuts == 0) = -1;  
    
    if nargin >= 3 % graphically show a tree        
        M = size(T,1);
        DG = sparse(1+T(:,1), 2:M+1, nMuts, M+1, M+1); 
        DG = DG(2:end, 2:end);
        bg = biograph(DG, ids, 'EdgeTextColor', [0.5 0 0.3], 'ArrowSize',3, ...
            'LayoutScale', 1.2, 'EdgeFontSize', 20); % dark red
        bg.ShowWeights = 'on';
        set(bg.Nodes,'Shape','circle');
        for t=1:M % go over all nodes in tree and add mutation string
            set(bg.Nodes(t), 'Description', mut_str{t});
        end                
        h = view(bg);
               
        set(gcf, 'Tag', '0');  
                
        h.ShowTextInNodes = 'label'; %shows the 'label' field as node label
        
        % Show a histogram for every need over the reads labels
        
        % time_stamp_data is an Mx2 matrix.  The first column is just the
        % node assignment of each read.  The second column is the label of
        % the read.
        
        % pie is M x nLabels matrix showing for each node how many reads 
        % come from each label        
        pie = T(:,2);  % default is all reads of the same one label

        % if we have a label for each sample we can add radio buttons
        % to show the induced tree for every single label        
        if nargin >= 4             
            btn = uibuttongroup('parent', gcf, 'Position',[0 0 1 0.05], 'Visible', 'Off');
            
            if pretty_max_colors > 1
                stamps = 0:pretty_max_colors;
            else
                stamps = [-1; unique(time_stamp_data(:,2))]; % list of labels 
                % (we prepended a -1 here to reserve a spot for the total
                % number of reads)
            end
                            
            % build the histogram for each node
            pie = zeros(M, length(stamps)); % if the max stamp was 3, then don't bother with 4,5,6,etc
            for u=length(stamps):-1:1  % for each label
                ix = time_stamp_data(:,2) == stamps(u);
                pie(:,u) = hist(time_stamp_data(ix,1), 1:M);
                nRun = sum(ix); % how many reads with that label total

                % create a radio button for that label
                cur_btn = uicontrol('Style','Radio','String',...
                    sprintf('label %d (%d reads)',stamps(u), nRun), ...
                    'pos',[10+150*(u-1) 5 145 15],'parent',btn,...
                    'HandleVisibility','off', 'tag', int2str(u-1), ...
                    'Visible', 'Off');%% Yi turned this off for figure generation reasons

            end
            
            % create the button that aggregates all labels
            pie(:,1) = sum(pie,2); % should be the same as T(:,2)
            set(cur_btn, 'String', sprintf('All (%d reads)', size(time_stamp_data, 1)));
            
            % Initialize some button-group (btn) properties. 
            set(btn,'SelectedObject',cur_btn);  % No selection
            set(btn,'SelectionChangeFcn',@resize_tree);
            set(btn,'Visible','Off'); %% Yi turned this off for figure generation reasons
            
            % encode in the figure the information required to draw the 
            % histogram of each node
            if ~isempty(sublabels)
                assert(M == size(sublabels,3), 'sublabels have the wrong size, as given to view_tree');
                for tt=1:M
                    h.Nodes(tt).Userdata.Distribution = sublabels(:,:,tt);
                end
            else % Joni's original logic for one-row format histogram
                for tt=1:M
                    h.Nodes(tt).Userdata.Distribution = pie(tt,2:end);
                end
            end
            
            
            bgInViewer.ShowTextInNodes = 'none';
            h.CustomNodeDrawFcn = @(node) dist_node_draw(node);
        end
        
        h.NodeAutoSize = 'off';
        resize_tree([], struct('NewValue', gcf)); 
        dolayout(h);
        h.scale = 0.65;
        if nargin >= 5
%            setAlwaysOnTop(figure(1),true); % not supported anymore.
            h.nodeCallbacks = @(node) show_node(node);
        end

    end
    
    % print the mutations in the sequences
    if nargin >= 2        
        for t=1:size(sequences,1)
            fprintf('%d:\t%s\n', t, char('0'+mut(t,:)));
        end
        max_mut = max(sum(mut>0));
        mut(mut == -2) = 5;        

        figure(2);
        subplot(3, 1, 1:2);
        view_read_alignment(mut); 
        ylim([0 max_mut+1]);
        title('summary of mutations');
        plot_class_lines( [1 (3-germline.frame_shift+1):3:size(mut,2) size(mut,2)+1], 0, 'k1', false);        
    end
    
    % function gets called when you click on a node
    function show_node(node)
        figure(1); 
        t = str2num(node.id);
        
        subplot(2, 1, 1);
        view_read_alignment(Z(time_stamp_data(:,1) == t, :), sequences(t,:));
        title('read noise')
        
        % show all mutations
        subplot(2, 1, 2);
        t_ = find(nodes == t);
        anc = X(t_,1:depth(t_)); 
        %anc = ancsestors(t, T);
        tmp = mut(anc, :); 
        tmp(1,:) = sequences(anc(1),:); % added July 09/2011
        if exist('germline', 'var')
            if ~isstruct(germline)
                germline = struct('seq', germline, 'frame_shift', 0);
            end
            %tmp(1,:) = 3+2*(germline.Sequence == 5); % so either 3 or 5
            map_('ACGTN')=1:5;
            tmp(1,:) = germline.seq; %map_(germline.Sequence);
        else 
            germline.frame_shift = 0;
            germline.seq = 5*ones(1, size(tmp,2));
        end        
        imagesc(tmp);        
        % frame_shift means:  how many bases were eaten from the beginning.
        plot_class_lines( [1 (3-germline.frame_shift+1):3:size(tmp,2) size(tmp,2)+1]);
        title('ancestral mutations up to germline'); 
        axis off
        
        return;
        tmp(tmp==5) = 0;
        tmp = tmp > 0;
        [G, L] = size(tmp);
        indexes = repmat((1:G)', 1, L);
        tmp = tmp.*indexes;
        for i=2:G
            tmp(i,:) = max(tmp(i,:), tmp(i-1, :));
        end
        subplot(3, 1, 2);
        bar(max(tmp));
        xlim([0 L]);        
        
        
    end

    % main function that control the node labels and node sizes
    function resize_tree(src, evt)        
        btn_chosen = evt.NewValue;
        tp = str2double(get(btn_chosen, 'tag'));

        new_size = pie(:,tp+1)';

        if ~isfield(evt, 'OldValue')
            old_size = inf;
        else
            old_chosen = evt.OldValue;
            tp_ = str2double(get(old_chosen, 'tag'));
            old_size = pie(:,tp_+1)';            
        end

        ix = (new_size == 0 & old_size ~= 0);
        set(h.Nodes(ix),'Label','','Size', [0.01 0.01]);

        ix = (new_size == 1 & old_size ~= 1);
        set(h.Nodes(ix),'Label','1','Size', [2 2]);
       
        ix = (new_size == 2 & old_size ~= 2);
        set(h.Nodes(ix),'Label','2', 'Size', [2 2]+1.5*sqrt(2));
        
        ix = find(new_size > 2);
        for t=ix
             set(h.Nodes(t), 'Label', int2str(new_size(t)));
             set(h.Nodes(t), 'Size', min([2 2]+1.5*sqrt(new_size(t)), [20 20]+log(new_size(t))*5));     
        end
        
        % Comment the section below to *not* show node id next to nodes.
        % We build mutation strings for the node e.g. 3A>T...                
        for t=2:length(h.Nodes) % t=1 corresponds to the germline node -- special description there
            % show '|'-delimited mutation strings on nodes vertically
            curr_node_mut_str = h.Nodes(t).Description; % Description contains full mutation string
            pipe_positions = find(curr_node_mut_str =='|');
            if length(pipe_positions) >= 2 % more than 2 delimiters
                ellipses_position = pipe_positions(1) - 1; % drop the final \n
                curr_node_mut_str = [curr_node_mut_str(1:ellipses_position) '...'];
            end
            curr_node_mut_str = regexprep(...
                curr_node_mut_str, '\|', '\n'); % leading part of mutstr is displayed
            
            set(h.Nodes(t),'Label', sprintf('%d [%d]\n%s',new_size(t),t, ...
                curr_node_mut_str));            
        end
        

        % replace -1 with 0
        edge_weights = get(h.Edges, 'weight');
        if length(edge_weights) == 1 % boundary -- singleton case
            if edge_weights == -1
                set(h.edges(1), 'weight', 0)
            end        
        else
            ix = find(cell2mat(edge_weights)== -1);
            set(h.edges(ix), 'weight', 0);
        end
    end     
end
