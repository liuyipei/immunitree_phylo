function [nuc2trim codon_positions123] = padded_nucseq_to_codon_position(...
  padded_query, padded_germline_inframe)

if nargin < 2
    unittest()
    codon_positions123 = [];
    fprintf('padded_nucseq_to_codon_position unit test passed!\n');
    return
end
% padded_nucseq_to_codon_position returns a vector containing 1,2,3
% indicating codon position. 0 refers to unknown.
% padded_query: sequence of 1xn nucleotides (in the form (1,2,3,4)), 
%       padded by skips in the form of numbers greater than 4.
%       This is the query sequence. 
% padded_VJ_inframe: sequence of 1nx nucleotides (in the form (1,2,3,4)),
%       padded by skips in the form of numbers greater than 4
%       This is the reference sequence (i.e. repertoire VJ)
% After padding, these sequences have the same length. The objective is to
% guess the frames of the bases in the query sequence. The first, second,
% and third bases of each codon are labeled 1,2,3 respectively. 0 is
% reserved for skips or unknowns.
% The return value is a sequence of numbers: 1,2,3 represent codon
% positions while 0 represents unknown.
%
% Example:
% padded_VJ_inframe         -ATCGA--TCG-TCAGT  (should be numeric)
% repVJ_indicator           01111100111011111
% repVJ_count               0123450067809ABCD
% repVJ_codon_position      01231200312031231
%                                   
% padded_query              AATCGAGGTCG-TC--T
% query_indicator           11111111111011001
% query_count               123456789AB0CD00E
% query_codon_position      12312312312031002
%                                  
% both_indicator            01111100111011001
% offsets3                  .11111..000.00..1
% offset               2                    
% codon_positions123        31231231231023001
%                           ^
%                           ^
% make inframe--> trim the '3' -- just 1 nucleotide

  assert(length(padded_query) == length(padded_germline_inframe))
  repVJ_indicator = padded_germline_inframe < 5;
  repVJ_count = cumsum(repVJ_indicator).*repVJ_indicator;
  repVJ_codon_position = mod(repVJ_count,3).*repVJ_indicator;

  query_indicator = padded_query < 5;
  query_count = cumsum(query_indicator).*query_indicator;
  query_codon_position = mod(query_count,3).*query_indicator;
  
  both_indicator = query_indicator & repVJ_indicator;
  
  % be sure to ignore positions that we can't say much about
  offsets3 = mod(3+repVJ_codon_position(both_indicator) - ...
      query_codon_position(both_indicator),3);  % short vector (neither seq can be padded)
  offset = mode(offsets3);
  
  nuc2trim = mod(3-offset, 3);
  
  if nargout < 2
      return
  end
  codon_positions123 = (1+mod(query_codon_position + offset+2,3)) .* query_indicator;   % making things 123, not 120
end

function unittest()
  map('ACGTNRacgtn-') = [1:5 5 1:5 5];
  padded_VJ_inframe = map('-ATCGA--TCG-TCAGT');
  padded_query =      map('AATCGAGGTCG-TC--T');
  [nuc_to_trim results] = ...
      padded_nucseq_to_codon_position(padded_query, padded_VJ_inframe);
  assert(nuc_to_trim == 1);
  assert(all(results == [3 1 2 3 1 2 3 1 2 3 1 0 2 3 0 0 1]));
end