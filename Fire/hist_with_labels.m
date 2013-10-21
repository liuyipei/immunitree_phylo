function [h bins] = hist_with_labels(x, labels, tl, yl, bins, normalize, group, cosmetic)
%function [h bins P H] = hist_with_labels(x, labels, tl, yl, bins, normalize, group, cosmetic)
    
    if ~exist('cosmetic', 'var'), cosmetic = []; end
    if ~exist('group', 'var'), group = 'grouped'; end
    if ~exist('normalize', 'var'), normalize = 0; end
    if ~exist('bins', 'var'), bins = 1; end
    if ~exist('labels', 'var') || isempty(labels), labels = ones(1,length(x)); end
    if ~isstruct(labels)
        l = labels;
        labels = [];
        labels.str = 'label';
        labels.num = max(l);
        labels.l = l;
        labels.names = strread(sprintf('%d\t',1:labels.num), '%s');
    end
    
    if isscalar(bins)    
        bins = min(x):bins:max(x);
    end

    if strcmp(group, 'grouped per patient')
        assert(isfield(labels, 'patients'));
        % For each patient show:
        % precentage of reads from each source with dist = X from germline
        figure;
        L = max(labels.patients);  assert(L == 10);
        h = zeros(length(bins), labels.num, L);
        ax = zeros(1,L);        
        for patient=1:L
            ax(patient) = subplot(L/2,2,patient);
            iz = labels.patients == patient;
            labels_ = labels;
            labels_.l = labels_.l(iz);
            h(:,:, patient) = hist_with_labels( x(iz), labels_, sprintf('%s (patient %d)',tl, patient), yl, bins, normalize, 'grouped');            
        end    
        linkaxes(ax, 'xy');
        return;
    end
        
    if strcmp(group, 'separate per patient')
        assert(isfield(labels, 'patients'));
        % For each patient show:
        % precentage of reads from each source with dist = X from germline
        figure;
        L = max(labels.patients);  assert(L == 10);
%        h = zeros(length(bins), labels.num, L);
        for patient=1:L
            iz = labels.patients == patient;
            labels_ = labels;
            labels_.l = labels_.l(iz);
            hist_with_labels( x(iz), labels_, sprintf('%s (patient %d)',tl, patient), yl, bins, normalize, 'cdfs', [L/2 2*labels.num (patient-1)*(labels.num)+1]);            
        end   
        return;
    end

    
    if exist('yl', 'var')
        if (normalize > 0)
            yl = ['% of ' yl];
        else
            yl = ['no. of ' yl];
        end
    else 
        yl = [];
    end
    
    if ~exist('tl', 'var')
        tl = '';
    end
    
    
    h = zeros(length(bins), labels.num);
    for l = 1:labels.num     % instances with label 0 are not shown in the plot
        ix = labels.l == l;
        h(:,l) = hist(x(ix & ~isnan(x)), bins)';        
    end

    if normalize == 1
        sumh = sum(h,1);
        h = 100*h ./ sumh(ones(size(h,1),1), :);
        tl_str = [' - ' labels.str ' normalized'];
    elseif normalize == 2
        h = 100*h / sum(h(:));
        tl_str = ' - total normalized';
    else 
        tl_str = '';
    end
    tl_str = [tl tl_str];
            
    if strcmp(group, 'cdfs') % with CDFs is is clear how we normalized
        tl_str = tl;
    end    

    if strcmp(group, 'cdfs')
        colors = 'rgcbm';
        if isempty(cosmetic)
            ax(l) = subplot(labels.num, 1, l);
        else
            ax(l) = subplot(cosmetic(1), cosmetic(2), cosmetic(3):(cosmetic(3)+labels.num-1));
        end
        for l = 1:labels.num % instances with label 0 are not shown in the plot
            stairs([0 bins], [0; cumsum(h(:,l))], colors(l), 'linewidth', 3);
            hold on;
        end
        title(tl_str);
        ylabel(yl);
        xlim([-0.5*(bins(2)-bins(1)) bins(ceil(4/5*length(bins)))]);
        ylim([0 100]);

    elseif strcmp(group, 'grouped')
        bar(bins, h, 'grouped');  
        title(tl_str);
        ylabel(yl);
        xlim([-0.5*(bins(2)-bins(1)) bins(ceil(4/5*length(bins)))]);
        ylim([0 max(h(:))]);

    elseif strcmp(group, 'separate')
        if cosmetic == 0, figure; end
        ax = zeros(1,labels.num);
        for l = 1:labels.num % instances with label 0 are not shown in the plot
            if isempty(cosmetic)
                ax(l) = subplot(labels.num, 1, l);
            else
                ax(l) = subplot(cosmetic(1), cosmetic(2), cosmetic(3) + l - 1);
            end
            bar(bins, h(:,l));            
            ylabel(yl);
            if l == 1
                title([tl_str ': ' labels.names{l}]);
            else
                title(labels.names{l});
            end
        end
        linkaxes(ax, 'xy');
        xlim([-0.5*(bins(2)-bins(1)) bins(ceil(4/5*length(bins)))]);
        ylim([0 max(h(:))]);
%         if isempty(cosmetic)
%             %figure;
%             [P H] = kstest_with_labels(x, labels, tl, 'aggregate');        
%         end
    else
        fprintf('error in ''group'' parameter\n');
    end    
    
end