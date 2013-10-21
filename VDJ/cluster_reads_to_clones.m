function [t, seqs, aligned, edit_str] = cluster_reads_to_clones(X, ...
    profileHMM_score, dup_map, nIterations, noise)
%function [t, seqs, aligned, edit_str] = cluster_reads_to_clones(X, profileHMM_score)
% cluster the reads in X to sequences.  Reads are given in integer form.
% output:
% t  cluster assignments.
% seqs sequences (in integers)
% aligned  - aligned version of each read (integers)
% edit_str - edit string for each read.

    if nargin == 0, [t, seqs, aligned, edit_str] = unittest; return; end
    alpha = 0.001;
    indel = 0.0015;
    if ~exist('noise', 'var'), 
        noise = 0.05;  % chance of mutation
    end
    nICM = 3;
    if ~exist('dup_map') || isempty(dup_map), dup_map = 1:length(X); end
    if ~exist('nIterations'), nIterations = 20; end
    N = length(dup_map);

    M = zeros(1,0); % counts how many reads in each cluster

    [~, ord] = sort(profileHMM_score(dup_map));
    ord = ord';
    assert(length(ord) == N);

    % Cache all alignments between reads/sequences
    nUnique = max(dup_map);
    cache = zeros(0,nUnique);
    cache_alignment = cell(0,nUnique);
    
    seqs = {};
    counts = {};    
    deletions = {};    
    trans = sparse(0,0);
    t = zeros(N,1);
    aligned = cell(N,1);
    dmap = cell(N,1);
    edit_str = cell(N,1);    
    
    for it = 1:nIterations        
        fprintf('iteration %d: ', it);
        moved = 0;
        cache_hits = 0;
        cache_tot = 0;
        trans(:) = 0;
        tic

        % Order reads according to likelihood of HMMs (best first) to encourage less clusters
        % for each read
        if it > 1, ord = randperm(N); end
        for i_=1:N
            i = ord(i_);
%            if i==16, keyboard; end
            if mod(i_, 1000) == 0, fprintf('%d ', i_); end

            % remove from stats
            old_ti = t(i); 
            if t(i) ~= 0
                M(t(i)) = M(t(i))-1;
                counts{t(i)} = counts{t(i)} - read_to_ind(aligned{i});
                deletions{t(i)} = deletions{t(i)} - (dmap{i}>0);
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
                    [cache(s,dup_map(i)) cache_alignment{s,dup_map(i)}] = ...
                        nw_alignment(seqs{s}, X{i}, indel, noise);
                end
            end
            ll_align = cache(:, dup_map(i))';
            ll = [log(M)+ll_align log(alpha)+profileHMM_score(dup_map(i))];
             % already take care of in above line +log(1-noise)*length(X{i})];

            if it > nIterations - nICM
                [~,k] = max(ll);
            else
                k = lmnrnd(ll, 0);
            end

            if k == length(ll)
                % new seq
                M = [M 0];
                seqs{k} = X{i};        
                counts{k} = zeros(5,length(X{i}));        
                deletions{k} = zeros(1,length(X{i})+1);        
                cache = [cache; +inf*ones(1,nUnique)];
                cache_alignment = [cache_alignment; cell(1,nUnique)];
                cache_alignment{k, dup_map(i)} = zeros(1,length(X{i}));
                trans(k,k) = 0;
            end   

            edit_str{i} = cache_alignment{k, dup_map(i)}; 
            %[i size(aligned) size(dmap) size(X,1) size(edit_str,1)]
            [aligned{i} dmap{i}] = apply_edit(X{i}, edit_str{i});

            % add to stats
            t(i) = k;
            M(t(i)) = M(t(i))+1;
            counts{t(i)} = counts{t(i)} + read_to_ind(aligned{i});
            deletions{t(i)} = deletions{t(i)} + (dmap{i}>0);
            if old_ti == 0 || ...
                    (t(i) ~= old_ti && M(t(i)) + M(old_ti) > 1)
                moved = moved + 1;
                if old_ti ~= 0
                    trans(old_ti, t(i)) = trans(old_ti, t(i)) + 1;
                end
            end

        end
        fprintf('...done. ');
        toc

        [seqs counts t M cache cache_alignment deletions trans] = clean_empty_seqs(seqs, ...
            counts, t, M, cache, cache_alignment, deletions, trans);

        
        % for each sequence
        old_seqs = seqs;
%        cache(:,end) = 0;  % all is clean
        for s = 1:length(M)
            % estimate consensus
            [~, seqs{s}] = max(counts{s}(1:4,:));
            % if sequence change, make dirty
            if ~isequal(old_seqs{s}, seqs{s})
%                cache(s,end) = +inf;
                cache(s,:) = +inf;
            end            
        end
        
        propose_indels = false;
        if propose_indels
        %  propose indels
        %  I already have a structure that stores insertions
        %  now how do I get one that stores deletions... done!
        %
            for s = 1:length(M)            
                thresh = M(s)/2;

                % Deletions
                fix = find(counts{s}(5, :) > thresh);
                if ~isempty(fix)
    %                cache(s,end) = -inf; % make dirty
                    cache(s,:) = +inf; % make dirty
                    for o=fliplr(fix)
                        % commit deletion in location o of seq s                    
                        fprintf('delete offset %d in sequence %d\n', o, s);
                        L = length(seqs{s});
                        del_fix = [1:o-1 o+1:L];
                        seqs{s} = seqs{s}(del_fix);
                        deletions{s} = zeros(1,L);
                        counts{s} = counts{s}(:, del_fix);
                        ix = find(t == s);
                        for i = ix'
                            dmap{i}(o+1) = dmap{i}(o) + dmap{i}(o+1) ...
                                + (aligned{i}(o) ~= 5);
                            dmap{i} = dmap{i}([del_fix L+1]);
                            aligned{i} = aligned{i}(del_fix);
                            deletions{s} = deletions{s} + (dmap{i}>0);
                            edit_str{i} = [];
                        end
                    end
                end

                % Insertions (more complicated?)
                fix = find(deletions{s} > thresh);
                if ~isempty(fix)
    %                cache(s,end) = -inf; % make dirty
                    cache(s,:) = +inf; % make dirty
                    fprintf('insert offsets [%s] in sequence %d\n', ...
                        num2str(fix), s);
                    for o=fliplr(fix)
                        seqs{s} = [seqs{s}(1:o-1) 5 seqs{s}(o:end)];
                    end
                    counts{s} = zeros(5,length(seqs{s}));
                    deletions{s} = zeros(1,length(seqs{s})+1);
                    ix = find(t == s);
                    for i = ix'
                        [~, edit_str{i}] = nw_alignment(seqs{s}, X{i}, indel, noise);
                        [aligned{i} dmap{i}] = apply_edit(X{i}, edit_str{i});
                        counts{s} = counts{s} + read_to_ind(aligned{i});
                        deletions{s} = deletions{s} + (dmap{i}>0);
                    end
                    [~, seqs{s}] = max(counts{s}(1:4,:));
                end

            end
        end
        

	% print some stats
        per_cache_hits =  100*cache_hits/cache_tot;
        fprintf('Total %d reads moved.  Total %d sequences %.1f%% cache hits.\n', moved, length(M), per_cache_hits);
	int16(M)
    int16(full(trans))
    end
end

function counts = read_to_ind(seq)
    L = length(seq);
    p = 0:5:(L-1)*5;
    counts = false(5,L);
    counts(seq+p) = 1;
end


function [seqs counts t M cache cache_alignment deletions trans] = ...
    clean_empty_seqs(seqs, counts, t, M, cache, cache_alignment, deletions, trans)
    % clean empty seqs
    ix = (M ~= 0);
    seqs = seqs(ix);
    counts = counts(ix);
    M = M(ix);
    cache = cache(ix, :);
    cache_alignment = cache_alignment(ix, :);
    deletions = deletions(ix);
    trans = trans(ix, ix);
    
    map = zeros(length(ix),1);
    map(ix) = 1:length(M);
    t = map(t);
end


% unit-test

function [t_, seqs_, Y, edt] = unittest()
%%
    % choose a V D J
    rep.V = fastaread('V.fa');
    rep.D = fastaread('D.fa');
    rep.J = fastaread('J.fa');
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
    for s=1:nSeq
        d = ceil(length(rep.D)*rand);
        D_seq = rep.D(d).Sequence; 
        [seqs{s} V N1 D N2 J ] = profileHMM_gen(V_seq, D_seq, J_seq); 
        seqs_nt{s} = map(seqs{s});
    end
    
    fprintf('Generating reads...\n');
    M = zeros(1,nSeq);    
    N = 100;
    noise = 0.02; indel = 0.003;
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
    [t_, seqs_, aligned, edt] = cluster_reads_to_clones(Y, profileHMM_score);
%    [t_' t]
%    confusion_matrix(t, t_)
%%
%     Y = cell(N,1);
%     map = 'ACGTN';
%     for i=1:N, Y{i} = map(X{i}); end
%     dist = seqpdist(Y,'Method','jukes-cantor',...
%                  'Indels','score',...
%                  'Alphabet', 'NT',...
%                  'PairwiseAlignment',true);

%            if isempty(alignments{k}) % cache hit but we don't have alignment
%                if k ~=old_ti 
%                    [~, edit_str{i}] = nw_alignment(seqs{k}, X{i});
%                    cache_hits = cache_hits - 1;
%                end
%            else
%                edit_str{i} = cache_alignment{k, dup_map(i)}; % cache miss
%            end


end

