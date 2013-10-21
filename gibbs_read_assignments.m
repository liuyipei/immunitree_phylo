function [t ll] = gibbs_read_assignments(tree, sequences, F, R, reads)

global greedy

logR = log(R); 

[unique_reads, ix, jx] = unique(reads, 'rows');
%unique_reads = (unique_reads-1)*4;

live_cells = find(tree(:,3) == F);
[unique_seqs, iy, jy] = unique(sequences(live_cells,:), 'rows');
unique_seqs = (unique_seqs-1)*5;
LLs = zeros(length(ix), length(iy));
for x=1:size(LLs,1)
    for y=1:size(LLs,2)
        iz = unique_seqs(y, :) + unique_reads(x, :);
        LLs(x,y) = sum(logR(iz));
    end
end
% uncomment if you want to sample from prior:
% LLs = zeros(length(ix), length(iy)); 

[t prob] = many_lmnrnd(LLs(jx,jy)',greedy);
ll = sum(LLs(sub2ind(size(LLs), jx, jy(t))));

t = live_cells(t);
ll = [ll ; -size(reads,1)*log(length(live_cells))];

end



