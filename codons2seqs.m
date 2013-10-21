function seqs = codons2seqs(codons, codon2nt)
if nargin<2
    global codon2nt
    %[~,codon2nt,~] = get_maps();
end
codon2nt_ = [codon2nt 5*ones(3,1)];

[N L] = size(codons);
seqs = zeros(N, 3*L, 'int16');
for i=1:N
    seqs(i,:) = reshape(codon2nt_(:,codons(i,:)), 1,[]);
end
end

function test()
%%
codons = [6 5 38 7; 30 47 10 24];
seqs = codons2seqs(codons);
assert(isequal(seqs2codons(seqs), codons));
fprintf('Test Passed.\n');

end