function [ll t] = slice_read_assignments(tree, sequences, t, F, R, reads, ...
    read_partition, annealing_prob)

% reminder:  columns of R sum to 1

% slice sampling

% need a method that will uniformly sample from the live cells

logR = log(R);

ix = read_partition.raw_from_uniq;
jx = read_partition.uniq_from_raw;
unique_reads = reads(ix,:);

% (over-stringent) correctness check against old way
%[unique_reads_old, ix_old, jx_old] = unique(reads, 'rows');
%assert(all(all(unique_reads_old == unique_reads)));
%assert(all(jx_old == jx));
%assert(all(ix_old == ix));


live_cells = find(tree(:,3) == F);

[X live_cells] = order_tree(tree, live_cells);
%clear X ord
map_live = zeros(max(live_cells),1);
map_live(live_cells) = (1:length(live_cells));
t = map_live(t);

% calculate current likelihood
[unique_seqs, iy, jy] = unique(sequences(live_cells,:), 'rows');
B = 5; 
unique_seqs = (unique_seqs-1)*B; % so that then we can "add"
LLs = zeros(length(ix), length(iy));


% for every read:
% calculate current likelihood
% calculate a threshold likelihood uniformly
% participate(read) = true.
N = size(reads,1);
participate = true(N,1);
proposed_likelihood = zeros(N,1);
interval = [ones(N,1), length(live_cells)*ones(N,1)];
for i=1:N
    x = jx(i); y = jy(t(i));
    if LLs(x,y) == 0
        LLs(x,y) = sum(logR(unique_seqs(y, :) + unique_reads(x, :)));
    end
    proposed_likelihood(i) = LLs(x,y);
end
thresh = log(rand(N,1))+proposed_likelihood;


if nargout == 2
    if true 
        % Bug fix -- revert SVN r91/92 -- better to diff this section
        % against SVN r31 by Joni. I still want to keep id in the same subclone.
        % -- Aug 28, 2012 (Yi)
        M = length(ix); % number of unique reads
        participate_uniq = true(M,1);
        
        interval_uniq = [ones(M,1), length(live_cells)*ones(M,1)];        
        t_uniq = t(ix);
        t_uniq_ = t_uniq;
        proposed_likelihood_uniq = sum(logR(unique_seqs(jy(t_uniq), :) + unique_reads),2);
        
        % scale by the number of duplicate sequences in thresh_uniq
        % scale by annealing prob
        thresh_uniq = log(max(rand(M,1), annealing_prob)) ./ read_partition.uniq_weights + proposed_likelihood_uniq(:);
        while any(participate_uniq)
        % for each column of the ancestors matrix (i.e. tree level).
        %   For each participating read, randomize a number between start_ix
        %    and end_ix of ancestor(l).  
            M_ = sum(participate_uniq);
            draw = rand(M_,1);
            t_uniq_(participate_uniq) = ...
                interval_uniq(participate_uniq,1) + ...
                floor(draw.*(interval_uniq(participate_uniq,2) - ...
                             interval_uniq(participate_uniq,1)+1));
        %   calculate new likelihood for each read
        %   if likelihood(read) > threshold(read): participate(read) = false.
            for x=find(participate_uniq)'                
                y = jy(t_uniq_(x));
                proposed_likelihood_uniq(x) = LLs(x,y);
                if LLs(x,y) == 0
                    LLs(x,y) = sum(logR(unique_seqs(y, :) + unique_reads(x, :)));
                    proposed_likelihood_uniq(x) = LLs(x,y);
                end
                participate_uniq(x) = proposed_likelihood_uniq(x) < thresh_uniq(x);
            end

            interval_uniq(participate_uniq & t_uniq_ < t_uniq,1) = ...
                       t_uniq_(participate_uniq & t_uniq_ < t_uniq);
            interval_uniq(participate_uniq & t_uniq_ > t_uniq,2) = ...
                       t_uniq_(participate_uniq & t_uniq_ > t_uniq);
        end
        t_ = t_uniq_(jx);

    else  % a different way of sampling
        X(X == 0) = -1;
        interval = make_intervals(X);
        t_ = t;
        for l=1:size(X,2)
            % for each column of the ancestors matrix (i.e. tree level).
            %   For each participating read, randomize a number between start_ix
            %    and end_ix of ancestor(l).
            participating = find(participate)';
            N_ = length(participating);
            if N_ == 0, break; end
            draw = rand(N_,1);
            anc = X(t(participating), l);
            ix = (anc == -1);
            t_(participating(ix)) = t(participating(ix));
            t_(participating(~ix)) = interval(anc(~ix),1) + floor(draw(~ix).*(interval(anc(~ix),2)-interval(anc(~ix),1)+1));
            %   calculate new likelihood for each read
            %   if likelihood(read) > threshold(read): participate(read) = false.
            for i=find(participate)'
                x = jx(i); y = jy(t_(i));
                proposed_likelihood(i) = LLs(x,y);
                if LLs(x,y) == 0
                    LLs(x,y) = sum(logR(unique_seqs(y, :) + unique_reads(x, :)));
                end
                proposed_likelihood(i) = LLs(x,y);
                participate(i) = proposed_likelihood(i)<thresh(i);
            end
        end
    end
    fprintf('%d reads moved to a different sequence\n', sum(jy(t_)~= jy(t)));
    t = live_cells(t_);
end

ll = [sum(proposed_likelihood) ; -N*log(length(live_cells))];
end

    
function intervals = make_intervals(X)

% for each column of the ancestors matrix (i.e. tree level).
%   do a diff
%   when it's not equal 0, put the node label.
%   find all non-zeros - that will give you the intervels.
%   add to table of ancs_id, start_ix, end_ix

intervals = zeros(size(X,2), 2);

junk_ = [ones(1,size(X,2)); diff(X)];
junk_(junk_~=0) = X(junk_~=0);
[I_,~, V_] = find(junk_);  % starting locations
ix = V_>0;
intervals(V_(ix),1) = I_(ix);

junk__ = flipud([ones(1,size(X,2)); diff(flipud(X))]);
junk__(junk__~=0) = X(junk__~=0);
[I__,~, V__] = find(junk__); % ending locations
ix = V__>0;
intervals(V__(ix),2) = I__(ix);


end


%     function res = ll_fun(x,y)
%         x = jx(x);
%         y = jy(y); %jy(map_live(y));
%         if LLs(x,y) == 0
%             iz = unique_seqs(y, :) + unique_reads(x, :);    
%             LLs(x,y) = sum(logR(iz));
%         end
%         res = LLs(x,y);
%     end


function test()
%%
t__ = best.t;
tic
[~,t__] = slice_read_assignments(best.tree, codons2seqs(best.sequences), t__, best.F, best.R, reads);  
%t__ = gibbs_read_assignments(best.tree, codons2seqs(best.sequences), best.F, best.R, reads);    
toc
sum(t__ ~= best.t)
end