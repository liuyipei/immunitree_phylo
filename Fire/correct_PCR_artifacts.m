function [data nSingle nNewSingle ix] = correct_PCR_artifacts(data, policy)

if nargin<2, policy = ''; end

nReads = length(data.reads);
M = length(data.subjects);
nSingle = 0;
nNewSingle = 0;

str = '  Correcting PCR artifacts...';
fprintf(['[' str ' '*ones(1,M-length(str)) ']\n[']);

% for each sample
for j=1:M
    fprintf('*');
    % for each PRIMER (FR1, FR2)
    for fr = 1:2

        % find unique reads
        ix = get_reads(data, j, [], fr, [], [0 0 0], 'ihmmune');
        [~,I J] = unique(data.reads(ix));
        N  = length(J);
        N_ = length(I);
        
        % construct 'counts':  For each unique read it has a row saying how
        % many reads from each replica
        tmp_counts = data.rep(ix)-'a'+1;
        counts = zeros(N_,max(tmp_counts));
        vals = 1:size(counts,2);
        for i=1:N_
            counts(i,:) = histc(tmp_counts(J==i), vals);
        end
        nSingle = nSingle + sum(sum(counts,2) == 1);

        if 0
            % all the unique reads that occur in a single replica   
            % leave only one occurance for each such read
            iy = find(sum(counts>0,2) == 1 & sum(counts,2) > 1);        
            for i=iy'
                iz = ix(J==i);
                data.rep(iz(2:end)) = 0;
            end
            nNewSingle = nNewSingle + length(iy);
        elseif 0
            % a more general approach:
            % take the minimum of those replicas where this read occurs
            % erase reads so that all replicas have the minimum                
            iy = find(sum(counts,2) > 1);
            for i=iy'
                iz = ix(J==i);
                [~,I_, J_] = unique(data.rep(iz));
                tmp = data.rep(iz(I_));
                data.rep(iz) = 0;
                data.rep(iz(I_)) = tmp;
            end        
            
            nNewSingle = nNewSingle + sum(sum(counts(iy,:)>0,2) == 1);
        else
            % erase reads from the replica has the maximum, until it has
            % the amount of the second-to-max plus one (if the max is
            % shared between 2 or more replicas, don't touch it
%            counts_ = counts;
            [maxes1, repmax] = max(counts,[],2);
            counts(sub2ind(size(counts), (1:N_)', repmax)) = 0;
            [maxes2, ~] = max(counts,[],2);
            iy = find(maxes1 > maxes2 + 1);
            for i=iy'
                newval = maxes2(i)+1;
                % all reads identical to i from the maximal replica
                iz = ix(J==i & tmp_counts == repmax(i));
                data.rep( iz( (newval+1):end) ) = 0;
            end

            nNewSingle = nNewSingle + sum((maxes1 > 1) & (maxes2 == 0));

        end
                
    end
end

% erase all records with rep = 0.

ix = data.rep ~= 0;

if strcmp(policy, 'removeSpam')
    ix = ix & (data.ihmmune(:,1) ~= 0);
end
    
flds = fieldnames(data); 
for k= 1:length(flds)
    if size(data.(flds{k}), 1) == nReads
        data.(flds{k}) = data.(flds{k})(ix,:);
    end
end

N_ = length(data.reads);



fprintf(']\n');
fprintf('\nTotal %d samples; %d reads.  %d (%.1f%%) PCR-duplicated reads removed.\n', M, nReads, nReads-N_, 100*(nReads-N_)/nReads);
fprintf('%.1f%% singletons before correction.\n%.1f%% singletons  after correction.\n', 100*nSingle/nReads, 100*(nSingle+nNewSingle)/N_);

end