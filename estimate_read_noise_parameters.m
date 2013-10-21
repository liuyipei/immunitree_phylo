function [R total_mismatches] = estimate_read_noise_parameters(sequences, reads, t, pseudo)

counts = zeros(5,5);

for u=1:size(reads,1)
    assert(t(u)>0);
    %data = mysub2ind(4, sequences(t(u),:), reads(u,:), 1);
    data = (sequences(t(u),:)-1)*5 + reads(u,:);
    % TODO:  remove this hist - it is computationally heavy
    counts = counts + reshape(histc(data, 1:numel(counts)), size(counts));        
end
total_mismatches = sum(counts(:))-sum(diag(counts));
R = generate_stochastic_matrix(pseudo + counts(1:5,:), 1, false);
end