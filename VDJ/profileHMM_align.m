function [score d alignments] = profileHMM_align(X, V, D, J, noise)

    % X is a cell array of reads
    N = length(X);
    nV = length(V);
    nD = length(D);
    nJ = length(J);

    % for now
    assert(nV == 1);
    assert(nJ == 1);
    
    params = profileHMM_get_params();
    if nargin==5, params.noise = noise; end
    
    profile_alignment = zeros(N, nD);
    a = struct('V', [], 'N1', [], 'D', [], 'N2', [], 'J', [], 'V_', [], 'D_', [], 'J_', [], 'eaten', [0 0 0 0], 'germline', []);
    a = a(ones(N,nD));
        
    for i=1:N 
        if mod(i, 100) == 0, fprintf('%d ', i); end
        if isempty(X{i}), continue; end
        for dix = 1:nD
            % align to profile HMMs
            if isempty(D{dix})
                [a(i, dix) profile_alignment(i, dix)] = ...
                  profileHMM_annotate_lambda(X{i}, V{1}, D{dix}, J{1}, params);
            else
                [a(i, dix) profile_alignment(i, dix)] = ...
                  profileHMM_annotate(X{i}, V{1}, D{dix}, J{1}, params);
            end
        end
    end
    if N>=100 fprintf('\n'); end
    
    [score d] = max(profile_alignment, [], 2);
    
    if nargout == 3
        alignments = a(sub2ind(size(a), 1:N, d));
    end
    
    % for sequences where there was no good alignment, output "d=0".
    d(isinf(score)) = 0;

end