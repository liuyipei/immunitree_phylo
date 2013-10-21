function Qs = get_Qs(mut_model)

%if nargin == 2, rate_class = 1; end

B = size(mut_model.NT,1);
L = length(mut_model.rate_class);
nClasses = size(mut_model.NT,3);

is_codon = isfield(mut_model, 'AA');
if is_codon
    assert(size(mut_model.AA,3) == nClasses);
    Qs = zeros(B^3, B^3, L);
else
    Qs = zeros(B, B, L);
end
   
is_decay = mut_model.decay ~= 0;

if is_decay || ~is_codon
    sum_NT = sum(mut_model.NT,1);
    NT_norm = mut_model.NT ./ sum_NT(ones(1,B),:,:);
end

if is_decay
    NT_delta = zeros(B, B, nClasses);
    for k=1:nClasses
        NT_delta(:,:,k) = NT_norm(:,:,k)^mut_model.decay;
    end   

    for l=1:L
        k = mut_model.rate_class(l);
        decayed_NT = real(NT_norm(:,:,k)*NT_delta(:,:,k)^l);
        decayed_NT(decayed_NT<0) = eps;
        if is_codon
            Qs(:,:,l) = make_codon_matrix(mut_model.AA(:,:,k), decayed_NT);
        else 
            sum_NT = sum(decayed_NT,1);
            Qs(:,:,l) = decayed_NT./sum_NT(ones(1,B),:);                        
        end
    end
else  % NO DECAY
    if is_codon
        for k=1:nClasses
            Qs(:,:,k) = make_codon_matrix(mut_model.AA(:,:,k), mut_model.NT(:,:,k));
        end   
    else
        Qs = NT_norm;
    end
    Qs = Qs(:,:,mut_model.rate_class);
end

Qs = reshape(Qs, size(Qs,1), []);
if abs(sum(Qs(:,1))-1)>1e-6, keyboard; end; 
if sum(Qs(:)<0) > 0, keyboard; end; 
end



function test()
%%
    AA = zeros(21, 21, 2);
    NT = zeros(4,4,2);

    AA(:,:,1) = 0.01+(1-21*0.01)*eye(21,21);
    NT(:,:,1) = 0.1+(1-4*0.1)*eye(4);

    AA(:,:,2) = 0.1+(1-21*0.1)*eye(21,21);
    NT(:,:,2) = 0.01+(1-4*0.01)*eye(4);

    Qs = get_Qs(AA,NT, [1 1 2]);
    len = numel(Qs)/3;
    assert(isequal(Qs(1:len),Qs((len+1):(2*len))));
    assert(~isequal(Qs(1:len),Qs((2*len+1):(3*len))));
    figure; imagesc(Qs);
end