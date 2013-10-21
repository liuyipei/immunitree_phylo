%% generate matrix of codons x mutation_offset x mut_value
function [SM I J K] = silent_map()
SM = false(65, 3, 4);
[codon2aa codon2nt nt2codon]= get_maps();

% This table summarizes what happens if letter number j of codon number i 
% becomes k (k=[1 2 3 4] for ACGT)
for i=1:64
    cod = codon2nt(:,i);
    aa = codon2aa(i);
    for j=1:3
        for k=1:4
            if cod(j) == k, continue; end  % it's already that letter
            cod_ = cod;            
            cod_(j) = k; % new codon
            aa_ = codon2aa(seqs2codons(cod_')); % new aa
            SM(i,j,k) = aa_ == aa;     % silent iff amino acid is the same
        end
    end
end



% assume SM (silent map) is in the workspace
% evaluates 3 64x64 tables.
I = false(64); % I(i,j)=1 iff codon(i)->codon(j) is a single-silent-mutation
J = false(64); % J(i,j)=1 iff codon(i)->codon(j) is a single-mutation 
K = false(64); % K(i,j)=1 iff codon(i)->codon(j) is a single-3rd-base-mutation 
for i=1:64
    for j=1:64
        x = codon2nt(:,i);
        y = codon2nt(:,j);
        ix = x ~= y;
        if sum(ix) == 1 % only consider single-mutations in codon
            J(i,j) = true;
            I(i,j) = SM(i, find(ix), y(find(ix)));
            K(i,j) = (find(ix) == 3);
        end
    end
end

% Hence to know if a mutation codon(i)->codon(j) is *not* silent we need to 
% have J(i,j) && ~I(i,j)
end