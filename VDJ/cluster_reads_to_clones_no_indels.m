function [t, seqs] = cluster_reads_to_clones_no_indels(X, profileHMM_score, dup_map, nIterations, noise)
%function [t, seqs] = cluster_reads_to_clones_no_indels(X, profileHMM_score, dup_map, nIterations, noise)
% cluster the reads in X to sequences.  Reads are given in integer form.
% output:
% t  cluster assignments.
% seqs sequences (in integers)

    if nargin == 0, [t, seqs] = unittest; return; end
    alpha = 0.1;
    
    if ~exist('noise', 'var'), 
        noise = 0.05;  % chance of mutation
    end
    log_noise = log(noise/3);
    log_match = log(1-noise);
    
     if ~exist('dup_map', 'var') || isempty(dup_map), 
         assert(size(X,1) == length(profileHMM_score));
         dup_map = 1:size(X,1);
     end
    
    N = length(dup_map);
    L = size(X,2);
    nICM = 2;
    X = double(X);
    
    if ~exist('nIterations'), nIterations = 20; end

    T = [1 1 1:0.5/(nIterations-nICM-3):1.5]; % temperature

    M = zeros(1,0); % counts how many reads in each cluster

    % Cache all alignments between reads/sequences
    nUnique = size(X,1);
    cache = zeros(0,nUnique);
    
    seqs = {};
    counts = {};    
    trans = sparse(0,0);
    t = zeros(N,1);
    
    % compute once the read indicators
    ind = zeros(5, L, N);
    for i=1:nUnique        
        ind(:,:,i) = read_to_ind(X(i,:));
    end
    
    for it = 1:nIterations        
        fprintf('iteration %d: ', it);
        moved = 0;
        cache_hits = 0;
        cache_tot = 0;
        trans(:) = 0;
        tic

        % Order reads according to likelihood of HMMs (best first) to encourage less clusters
        % for each read
        ord = randperm(N);
        for i_=1:N
            i = ord(i_);
            if mod(i_, 1000) == 0, fprintf('%d ', i_); end

            % remove from stats
            old_ti = t(i); 
            if t(i) ~= 0
                M(t(i)) = M(t(i))-1;
                counts{t(i)} = counts{t(i)} - ind(:,:,dup_map(i));
                t(i) = 0;
            end        

            % for each existing sequence
            for s = 1:length(M)
                if M(s) == 0, cache(s,dup_map(i)) = -inf; continue; end
                % if seq s is clean, get score from cache.
                cache_tot = cache_tot + 1;
                if cache(s,dup_map(i)) < 0 % cache hit
                    cache_hits = cache_hits+1;
                else
                    x_i = X(dup_map(i),:);
                    nMatch = sum(seqs{s}==x_i);
                    nGaps = sum(x_i == 5);
                    cache(s,dup_map(i)) = log_noise*(L-nMatch-nGaps) + log_match*nMatch;
                end
            end
            ll_align = cache(:, dup_map(i))';
            ll = [log(M)+ll_align log(alpha)+profileHMM_score(dup_map(i))];
             % already taken care of in above line +log(1-noise)*length(X{i})];

            if it > nIterations - nICM
                [~,k] = max(ll);
            else
                k = lmnrnd(T(it)*ll, 0);
            end

            if k == length(ll)
                % new seq
                M = [M 0];
                seqs{k} = X(dup_map(i),:);
                counts{k} = zeros(5,L);        
                cache = [cache; +inf*ones(1,nUnique)];
                trans(k,k) = 0;
            end   

            % add to stats
            t(i) = k;
            M(t(i)) = M(t(i))+1;
            counts{t(i)} = counts{t(i)} + ind(:,:,dup_map(i));
            if old_ti == 0 || ...
                    (t(i) ~= old_ti && M(t(i)) + M(old_ti) > 1)
                moved = moved + 1;
                if old_ti ~= 0
                    trans(old_ti, t(i)) = trans(old_ti, t(i)) + 1;
                end
            end

        end
        fprintf('...done. ');

        [seqs counts t M cache trans] = clean_empty_seqs(seqs, ...
            counts, t, M, cache, trans);

        
        % for each sequence
        old_seqs = seqs;
        for s = 1:length(M)
            if it > nIterations - nICM
                % estimate consensus            
                [~, seqs{s}] = max(counts{s}(1:4,:));
            else
                ll = counts{s}*(log_match-log_noise);
                seqs{s} = many_lmnrnd(T(it)*ll(1:4,:), 0);
            end
            % if sequence change, make dirty
            if ~isequal(old_seqs{s}, seqs{s})
                cache(s,:) = +inf;
            end            
        end
                          
        % print some stats
        per_cache_hits =  100*cache_hits/cache_tot;
        fprintf('Total %d reads moved.  Total %d clones, largest size %d.  %.1f%% cache hits.\n', moved, length(M), max(M), per_cache_hits);
    end
end

function counts = read_to_ind(seq)
    L = length(seq);
    p = 0:5:(L-1)*5;
    counts = false(5,L);
    counts(seq+p) = 1;
end

function total = multi_read_to_counts(X, B)
    total = histc(X, 1:B, 1);

    X = double(X);    
    total_ = zeros(B, size(X,2));
    ix = 0:B:(size(X,2)-1)*B;
    for k=1:size(X,1)
        iz = ix+X(k,:);
        total_(iz) = total_(iz)+1;
    end
end

function [seqs counts t M cache trans] = ...
    clean_empty_seqs(seqs, counts, t, M, cache, trans)
    % clean empty seqs
    ix = (M ~= 0);
    seqs = seqs(ix);
    counts = counts(ix);
    M = M(ix);
    cache = cache(ix, :);
    trans = trans(ix, ix);
    
    map = zeros(1,length(ix));
    map(ix) = 1:length(M);
    t = map(t);
end


% unit-test

function [t_, seqs_, Y, edt] = unittest()
%%
    % choose a V D J
    rep.V = fastaread('rep/IGH/V.fa');
    rep.D = fastaread('rep/IGH/D.fa');
    rep.J = fastaread('rep/IGH/J.fa');
    % choose V, D, J.
    v = ceil(length(rep.V)*rand);
    j = ceil(length(rep.J)*rand);
    V_seq = rep.V(v).Sequence; 
    J_seq = rep.J(j).Sequence;

%%    
    map('ACGTNacgtn') = [1:5 1:5];
    nSeq = 10;
    seqs = cell(nSeq,1);
    fprintf('Generating sequences...\n');
    d = ceil(length(rep.D)*rand);
    for s=1:nSeq        
        D_seq = rep.D(d).Sequence; 
        [seqs{s} V N1 D N2 J] = profileHMM_gen(V_seq, D_seq, J_seq);
        seqs_nt{s} = map(seqs{s});
    end
    
    fprintf('Generating reads...\n');
    M = zeros(1,nSeq);
    N = 100;
    noise = 0.02; indel = 0; %0.003;
    [X t] = gen_reads(seqs, N, noise, indel);

%%   
    fprintf('Trimming reads to CDR region...\n');
    [Y] = trim_reads(X, V_seq, J_seq);
    
%%    
    fprintf('Mapping the reads to profile HMMs...\n');
    tic;
    [profileHMM_score profileHMM_id] = profileHMM_align(Y, {V_seq}, {rep.D.Sequence}, {J_seq});
    toc;

%%
    ix = find(t==2 | t==9);
    %Y
    %profileHMM_score    
    %cell2mat(Y(ix))
    profileHMM_score(ix)
    [t_ seqs_] = cluster_reads_to_clones_no_indels(cell2mat(Y(ix)), profileHMM_score(ix));

end