function Z = plot_patient_repertoire(big_hist,Vs, Js, ttl)
    
    if ~exist('ttl', 'var'), ttl = []; end;
    if ~exist('Vs', 'var'), Vs = []; end;
    if ~exist('Js', 'var'), Js = []; end;
        
    if size(big_hist,3) == 1
        imagesc(big_hist');
        if ~isempty(Vs)
            h = text(1:length(Vs), 0.5+1*(-1+mod(1:length(Vs),2)), ...
                cellfun(@(x) x(4:end), Vs, 'uniformoutput', false));
            set(h, 'rotation', 90);
            set(h, 'fontsize', 8);
        end

        if ~isempty(Js)
            h = text(2*(-1+1*mod(1:length(Js),2)), 0.25+(1:length(Js)), cellfun(@(x) x(4:end), Js, 'uniformoutput', false));
%            h = text(zeros(1,length(Js)), 0.25+(1:length(Js)), cellfun(@(x) x(4:end), Js, 'uniformoutput', false));
            set(h, 'HorizontalAlignment', 'right')
        end
    %    set(gca, 'position', [0.11 0.11 0.775 0.515]);
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])

        if ~isempty(ttl)
            h = text(size(big_hist,1)+9,3, ttl);
            set(h, 'HorizontalAlignment', 'left')
        end
    %    xlabel('all patients together');
        colorbar
        colormap hot
        Z = big_hist;
        
    else
        % we assume that big_hist is 3 dimensional and each sample's
        % histogram is serialized in the first dimension.
        Z = [];
        for k=1:size(big_hist,3)
            Z = [Z ; reshape(big_hist(:, :, k), length(Js), [])];
        end
        if ~isempty(ttl)
            maxZ = max(Z(:));
            %imagesc(Z, [0 0.1]);
            imagesc(Z, [0 min(maxZ,600)]);
            
            colormap hot;
            plot_class_lines(1:length(Vs):size(Z,2)+1, false, 'm', true);
            plot_class_lines(1:length(Js):size(Z,1)+1, true, 'm', true);
            title(ttl) 
            pos = get(gcf, 'position');
            pos(3:4) = [1024 768];
            set(gcf, 'position', pos);

            if iscell(Vs)
                Vs = repmat(Vs, 1, 4);
                h = text(1:length(Vs), 0.5+1*(-1+mod(1:length(Vs),2)), ...
                    cellfun(@(x) x(4:end), Vs, 'uniformoutput', false));
                set(h, 'rotation', 90);
                set(h, 'fontsize', 8);

                set(gca, 'xtick', 1:3:(57*4))
                set(gca, 'xticklabel', 1:3:57)
                set(gca, 'XAxisLocation', 'top');

                
%                 h = text(1:length(Vs), 65.5+1*(-1+mod(1:length(Vs),2)), ...
%                     cellfun(@(x) x(4:end), Vs, 'uniformoutput', false));
%                 set(h, 'rotation', 90);
%                 set(h, 'fontsize', 8);
            end

        end        
    end
end
    