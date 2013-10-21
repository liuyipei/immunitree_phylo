function [codon2aa codon2nt nt2codon]= get_maps()
    dict = 'ACGT';
    [NT1 NT2 NT3] = ind2sub([4 4 4], 1:64);
    codon2nt = int16([NT1; NT2; NT3]);
    nt2codon = int16(1:64);
    nt2codon  = reshape(nt2codon, 4, 4, 4);
    codon2aa = int16(aa2int(nt2aa(dict(codon2nt'))));
    codon2aa(codon2aa == 24) = 21;
%    codon2aa(65) = 22;  % unknown?
end
