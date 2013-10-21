function [Y action] = trim_reads(X, V, J)
% Obselete - don't use.  Use trim_and_fix_reads instead.
%function Y = trim_reads(X, V, J)
% Given a set of reads X, and a prefix V, and a suffix J, returns a new seq
% of reads that starts exactly where the V starts, and end exactly where
% the J ends.

    N = length(X);
    st = zeros(2, N);
    en = zeros(2, N);
    Y = cell(N,1);
    
% preperation    
    eps = 0.001;
    indel = 0.003;
    sigma = 0.02;
    
    p_match_to_unobserved_seq = log(0.25);
    p_match_to_unobserved_read = 0;
    p = (1-4*sigma/3)*eye(5);
    p = p + sigma/3;
    p = log(p);
    p(end,:) = p_match_to_unobserved_seq;
    p(:, end) = p_match_to_unobserved_read;
    q = log(0.25);
    penalty = -(log(indel)+q);
    corr = 6;
    map('ACGTNacgtn-') = [1:5 1:5 5];
    V = map(V);
    revJ = fliplr(map(J));
    
    for i=1:N

        if mod(i,100) == 0, fprintf('%d ', i); end
        if X{i}(1) > 5, X{i} = map(X{i}); end
    
        if ~isempty(V)
            % a better thing to do would be:  align V to X{i} using SW.
            [s al st(:,i)] = swalign(V, X{i}, 'alphabet', 'nt', 'ScoringMatrix', p+corr, ...
                'GapOpen', penalty, 'ExtendGap', penalty-eps-corr);
        end
        
        if ~isempty(J)
            % align J to X{i}.
            [s a en(:,i)] = swalign(revJ, fliplr(X{i}), 'alphabet', 'nt', 'ScoringMatrix', p+corr, ...
                'GapOpen', penalty, 'ExtendGap', penalty-eps-corr);
        end
    
    end
    fprintf('...done.\n');    
    action = [st(2,:)-st(1,:); en(2,:)-en(1,:)]';
    for i=1:N
        % Find what matches to first letter of V
        if action(i,1) >= 0 % delete prefix of read
            Y{i} = X{i}(action(i,1)+1:end);
        else             % add insertions to beginning of read
            Y{i} = [5*ones(1,-action(i,1)) X{i}];
        end
        
        % Find what matches to last letter of J.
        if action(i,2) >= 0 % delete suffix of read
            Y{i} = Y{i}(1:end-action(i,2));
        else              % add insertions to end of read
            Y{i} = [Y{i} 5*ones(1,-action(i,2))];
        end                
    end
end



function unittest()
%%    
    rep.V = fastaread('V.fa');
    rep.J = fastaread('J.fa');
    v = 6; 
    j = 4;
    V_seq = rep.V(v).Sequence; 
    J_seq = rep.J(j).Sequence;

    Z = fastaread(sprintf('/afs/cs/u/joni/scratch/data/Uri/bins/bin_%d_%d.fasta', v, j));
    X = {Z.Sequence};
%%    
    [Y] = trim_reads(X, V_seq, J_seq);

end


