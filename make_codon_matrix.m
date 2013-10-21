function Q = make_codon_matrix(AA, NT)
    global codon2aa codon2nt
    Q = AA(codon2aa,codon2aa);
    for n=1:3
        Q = Q.*NT(codon2nt(n,:), codon2nt(n,:));
    end
    sumQ = sum(Q,1);  % Q is (child,parent)
    Q = Q./sumQ(ones(1,64),:);    
end