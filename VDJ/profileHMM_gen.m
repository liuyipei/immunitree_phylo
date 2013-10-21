function [seq V N1 D N2 J]  = profileHMM_gen(v, d, j, params)
    
    if nargin < 4
        params = profileHMM_get_params();
    end
   
    % choose end-deletions
    for i=1:4
       p(i) = lmnrnd(log(params.P_del(i).pdf), 0); 
    end

    % choose N-additions
    N = []; n = [];
    for i=1:2
       n(i) = lmnrnd(log(params.N_add(i).pdf), 0); 
       N{i} = super_lmnrnd(log(params.N_emit), n(i));
    end

    % make base sequence:
%    map = [];
%    map('ACGTNacgtn') = [1:5 1:5];

    % cut necessary parts from V
    assert(length(v) >= p(1));
    V = v(1:(end-p(1)));

    % cut necessary parts from D
    if p(2) + p(3) > length(d)
        D = [];
    else
        D = d( (p(1)+1):(end-p(2)) );
    end

    % cut necessary parts from J
    assert(length(j) >= p(4));
    J = j((p(4)+1):end);

    % generate seq
    map = 'ACGTN';
    N1 = map(N{1}); 
    N2 = map(N{2});
    
    ix = find(V == 'N' | V == 'n'); V(ix) = map(ceil(4*rand(1,length(ix))));
    ix = find(D == 'N' | D == 'n'); D(ix) = map(ceil(4*rand(1,length(ix))));
    ix = find(J == 'N' | J == 'n'); J(ix) = map(ceil(4*rand(1,length(ix))));

    seq = [V N1 D N2 J];

%     % fix seq for the case the germline has N nucleotides in it
%     ix = find(seq == 5);
%     seq(ix) = ceil(4*rand(1,length(ix)));
end