function [stop_codon_indices] = find_stop_codons_given_unpadded_inframe(...
  aStopCodonPerRow, inputNucleotides)
    stop_codon_indices = []; % default return value is none are found
    numTrios = floor(length(inputNucleotides) / 3);
    for i = 1:numTrios
        seq_indices = (i*3) - [2 1 0];
        for j = 1:size(aStopCodonPerRow,1)
            curr_input = inputNucleotides(seq_indices);
            curr_stop = aStopCodonPerRow(j,:);
            if all(curr_stop== curr_input); % match all 3 nucs
                stop_codon_indices = seq_indices;
                return
            end
        end
    end
    
end
