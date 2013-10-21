% convert seqs to codons, changing each triplet to a number between 1 to 64
function [codons AAs] = seqs2codons(seqs)
    [M L] = size(seqs);
%    if mod(size(L),3) ~= 0, L = L-mod(size(L),3); end
    assert(mod(L,3) == 0);
    seqs = reshape(seqs', 3,[]);
    codons = mysub2ind(4, seqs(1,:), seqs(2,:), seqs(3,:));
    ix = (seqs(1,:) == 5) | (seqs(2,:)==5) | (seqs(3,:) == 5);
    codons(ix) = 65;
    codons = reshape(codons, L/3, M)';
    if nargout > 1
        global codon2aa; %codon2aa = get_maps();        
        codon2aa_ = [codon2aa; 22];
        AAs = codon2aa_(codons);
    end
end



function test()
%%
seqs = [3      3      2      4      3      3      1      3      4;
      3      2      2      4      1      3      1      3      4;
      3      2      2      4      1      3      1      3      4;
      3      3      4      2      4      1      3      1      4;
      3      2      2      4      1      3      1      3      4;
      3      3      2      4      3      3      1      3      4;
      3      2      2      4      1      3      1      3      4;
      3      2      2      4      1      3      5      3      4;
      3      2      2      4      1      3      1      3      4;
      3      2      2      4      1      3      1      3      4];
[codons AAs] = seqs2codons(seqs)


end

