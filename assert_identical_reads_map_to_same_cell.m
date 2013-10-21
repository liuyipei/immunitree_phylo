
function assert_identical_reads_map_to_same_cell(raw_from_uniq, uniq_from_raw, t_)
    %[~, raw_from_uniq, uniq_from_raw] = unique(reads, 'rows');

    uniq_reads_t_ = zeros(1,length(raw_from_uniq));
    for i = 1:length(uniq_from_raw)
        x = uniq_from_raw(i);
        if uniq_reads_t_(x) == 0
            uniq_reads_t_(x) = t_(i);
        else
            assert(uniq_reads_t_(x) == t_(i), 'Error: identical reads were mapped to different cells!')
        end
    end
end