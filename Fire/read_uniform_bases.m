function best = read_uniform_bases(best, Vt_pos_diverse, Vt_mode)    
    best.sequences_diverse = best.sequences;
    best.sequences = zeros(size(best.sequences_diverse, 1), length(Vt_mode));
    best.sequences(:, Vt_pos_diverse) = best.sequences_diverse;
    best.sequences(:,~Vt_pos_diverse) = repmat(Vt_mode(~Vt_pos_diverse), ...
        size(best.sequences_diverse, 1) ,1);
end
