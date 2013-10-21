% compares the sequences in Y (one line per sequence) to the sequence
% 'seq'.  Return 'Y' which shows the mutations only, and total, which just
% counts how many mutation of each letter were in each column.
% If there are no output parameters, shows it in a graph.
function [Y total] = view_read_alignment(Y, seq)
if nargin<2
    seq = 5*ones(1,size(Y,2));    
%    seq = mode(Y);
    
end

T = [0; ones(size(Y,1),1)];
[~,Y] = annotate_mutations_on_tree(T, [seq; Y], true);

% Y(Y==-2) = 5; 
% -2 means that there was no mutation with respec to the germline
Y(Y==-2) = 7; 
Y = Y(2:end,:);

total = zeros(5, size(Y,2));
for k=1:5
    total(k,:) = sum(Y == k, 1);
end
if nargout == 0
    total_ = [sort(total(1:4, :)); total(5,:)];
    bar_handle = bar(total_', 1, 'stacked');
    if ~isempty(Y)
        xlim([0 size(Y, 2)]);
        ylim([0 min(size(Y, 1), max(sum(total(1:4,:)))+1)]);
        ylabel(sprintf('total %d', size(Y,1)));
    end
end

end

