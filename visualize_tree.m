function h = visualize_tree(T, seqs, t, X, nLabels, germline, labels, sublabels)

if ~exist('nLabels', 'var'), nLabels = 1; end
if ~exist('germline', 'var'), germline.seq = seqs(1,:); germline.frame_shift = 0; end


assert(~isempty(t));

if size(T,1) >1
    if size(t,2) == 2        
        if ~exist('X', 'var')        
            h = view_tree(T, seqs,nLabels, t);            
        else
            h = view_tree(T, seqs,nLabels, t, X, germline);            
        end
    else
        if ~exist('X', 'var')
            view_tree(T, seqs,1);
        elseif ~exist('labels', 'var')  % error rate representation
            labels = zeros(size(t(:)));
            dict = (0:size(seqs,2))';
            dict(5:10) = 4; % 4 or more errors
            dict(11:20) = 5; % 10 or more errors
            dict(21:end) = 6; % 20 or more errors
            pretty_max_colors = 6;
            for i=1:size(X,1)                
                [nw_seq1, nw_seq2, ~] = prettyalign2seqs(X(i,:), seqs(t(i),:));
                labels(i) = sum(nw_seq1 ~= nw_seq2);                
            end
            labels = dict(labels + 1) + 1;
            
            h = view_tree(T, seqs, pretty_max_colors, [t labels], X, germline);
        else % for example, if I wish to color the histograms based on time
            if exist('sublabels', 'var') % multi visit, multi isotypes, leading to two-layered pie-charts
                pretty_max_colors = 0;
                h = view_tree(T, seqs, pretty_max_colors, [t labels], X, germline, sublabels);
            else
                % single-layered pie-charts, basically just a single integer index associated
                % with each read , just like the error-rate representation                
                pretty_max_colors = max(6, max(labels));
                h = view_tree(T, seqs, pretty_max_colors, [t labels], X, germline);
            end
        end
    end
else 
    h = [];
end
end
function [seq1, seq2, aligned] = prettyalign2seqs(seq1, seq2)
    % convert numbers into letters, do alignment, and convert back into numbers
    map_from_int = 'ACGT-';
    map_to_int = [];
    map_to_int('ACGTNRacgtn-') = [1:5 5 1:5 5];
    seq1_ = map_from_int(non5(seq1));
    seq2_ = map_from_int(non5(seq2));
    [~, aligned] = nwalign(seq1_, seq2_);
    seq1 = map_to_int(aligned(1,:));
    seq2 = map_to_int(aligned(3,:)); % second row is colons, pipes, spaces
end

function non5 = non5(x)
    non5 = x(x~=5);
end