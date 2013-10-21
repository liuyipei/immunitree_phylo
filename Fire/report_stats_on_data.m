%% Report on the reads
fprintf('\n\t\ttotal\tFR1\tFR2\n---------------------------------------\n');
fprintf('total\t\t%6d\t%6d\t%6d\n', length(data.reads), sum(data.FR == 1), sum(data.FR == 2));
fprintf('in-frame\t%6d\t%6d\t%6d\n', sum(data.RF == 10), sum(data.FR == 1 & data.RF == 10), sum(data.FR == 2 & data.RF == 10));
fprintf('out code 20\t%6d\t%6d\t%6d\n', sum(data.RF == 20), sum(data.FR == 1 & data.RF == 20), sum(data.FR == 2 & data.RF == 20));
fprintf('out code 30\t%6d\t%6d\t%6d\n', sum(data.RF == 30), sum(data.FR == 1 & data.RF == 30), sum(data.FR == 2 & data.RF == 30));
fprintf('out code 40\t%6d\t%6d\t%6d\n', sum(data.RF == 40), sum(data.FR == 1 & data.RF == 40), sum(data.FR == 2 & data.RF == 40));

hist(data.RF+data.FR, 1:50);
%% How many reads have gaps
gap_counter = 0;
gap_hist = zeros(1,10);
for i=1:length(data.reads)    
    if mod(i,10000) == 0, fprintf('%d...', i); end
    gap_read = sum(data.reads{i} == 'N');
    %assert(gap_read == sum(nt2int(data.reads{i}) > 4)); % didn't fail
    if gap_read > 0
        gap_counter = gap_counter+1;
        gap_hist(gap_read) = gap_hist(gap_read) + 1;
    end    
end
fprintf('\n');
gap_counter 
gap_hist

%% Analyze multiple V combinations
N= length(data.reads);
vdj_hist = hist(data.ihmmune(:,[1 4 5 6 7]), 0:max(data.ihmmune(:)));
vdj_hist(1,2:end) = 0; 
subplot(4,1,1:2);
bar(vdj_hist, 'stacked');
title('V alignment popularity');
fprintf('%.1f%% of all reads do not have V alignment\n', 100*sum(data.ihmmune(:,1) == 0)/N);

% Analyze how many reads do not have a D alignment.
subplot(4,1,3);
hist(data.ihmmune(:,2), 0:max(data.ihmmune(:,2)));
title('D alignment popularity');
fprintf('%.1f%% of all reads do not have D alignment\n', 100*sum(data.ihmmune(:,2) == 0)/N);

% Analyze how many reads do not have a J alignment.
subplot(4,1,4);
hist(data.ihmmune(:,3), 0:max(data.ihmmune(:,3)));
title('J alignment popularity');
fprintf('%.1f%% of all reads do not have J alignment\n', 100*sum(data.ihmmune(:,3) == 0)/N);


%% Analyze replicates
figure; hist(data.rep-'a'+double(10*data.FR), 10:24);
title('breakdown of reads into replicates (left - FR1, right - FR2)');

%% uniqueness?
reads = data.reads;
ix = find(data.FR == 1);
for i=ix'    
    reads{i} = reads{i}(64:end);
end;
Z = cellfun(@length, reads);
hist(Z, 1:400);
title('length histogram of reads (FR1 reads trimmed by 63b)');
%%
%im = 1:4:40;  im = [im+1; im+2; im+3]; im = im(:)
im = 1:40;
M = length(im);
stats = zeros(M,8);

for j=1:M
    j
    ix = get_reads(data, im(j), [], [], [], [0 0 0], 'ihmmune');
    [~,I J] = unique(data.reads(ix));

    N  = length(J);
    N_ = length(I);
    counts = zeros(N_,6);
    tmp_counts = 10*double(data.FR(ix))+data.rep(ix)-'a';
    vals = [10 11 20 21 22 23];    
    for i=1:N_
        counts(i,:) = histc(tmp_counts(J==i), vals);
    end
    % count no. reads and no. distinct reads
    stats(j, 1:2) = [N N_];
    
    % no. of singletons
    stats(j,3) = sum(sum(counts,2) == 1);
    
    % no. of non-singletons but in a single replica
    iy = sum(counts>0,2) == 1 & sum(counts,2) > 1;
    stats(j,4) = sum(iy);
    stats(j,5) = mean(sum(counts(iy,:),2));  % mean number of duplicates

    % no. of non-singletons and in more than one replica
    iy = sum(counts>0,2) > 1;    
    stats(j,6) = sum(iy);    
    stats(j,7) = mean(sum(counts(iy,:),2)); % mean no. of duplicates
    iy = counts(iy,:);
    stats(j,8) = mean(iy(iy>0)); % mean no. of duplicates per replica
    
end
%%
fprintf('\ntotal\t     distinct\t   singletons\t    non-single 1 replica\tnon-single >1 replica\n');%-----------------------------------------------------------------------------\n');
fprintf('       \t             \t             \t        (avg per replica)\t(avg total, avg per replica)\n-----------------------------------------------------------------------------------------------\n');
for j=1:M
    fprintf('%6d\t\t%2.1f%%\t\t%2.1f%%\t\t%2.1f%% (%.1f)\t\t%2.1f%% (%.1f, %.1f)\n', stats(j,1), 100*stats(j,2)/stats(j,1), 100*stats(j,3)/stats(j,1), 100*stats(j,4)*stats(j,5)/stats(j,1), stats(j,5), 100*stats(j,6)*stats(j,7)/stats(j,1), stats(j,7), stats(j,8));
    if mod(j,4) == 0, fprintf('\n'); end
end


%%  DONE


%% Reads only appearing once

N  = length(J);
N_ = length(I);
fprintf('\nThere are %d unique reads (%.1f%% of the total %d reads).\n', N_, 100*N_/N, N);
fprintf('Of those:\n');
fprintf('\t%6d (%.1f%%) appear only once (%.1f%% of all reads).\n', sum(sum(counts,2) == 1),100*sum(sum(counts,2) == 1)/N_ , 100*sum(sum(counts,2) == 1)/N );
%fprintf('%d of them (%.1f%%) appear only once (%.1f%% of all reads).\n', sum(sum(counts,2) == 1),100*sum(sum(counts,2) == 1)/N_ , 100*sum(sum(counts,2) == 1)/N );

ix = sum(counts>0,2) == 1 & sum(counts,2) >= 1;
h = sum(counts(ix,:),2);
h = hist(h, 1:max(h));
fprintf('\t%6d (%.1f%%) appear more than once, but on the same replicate (%.1f%% of all reads).\n', sum(h)-h(1),100*(sum(h)-h(1))/sum(h), 100*h(2:end)*(2:length(h))'/N);
subplot(3,1,1);
h = sum(counts(ix,:),2);
h = hist(h, 1:max(h));
bar(h/sum(h));
title('reads in only one replicate');

fprintf('\t\tOf those:\n');
ix = sum(counts>0,2) == 1 & sum(counts(:,1:2),2) >= 1;
h1 = sum(counts(ix,:),2);
h1 = hist(h1, 1:max(h1));
fprintf('\t\t%5d reads (%.1f%%) in FR1, %.1f duplicates on average\n', (sum(h1)-h1(1)), 100*(sum(h1)-h1(1))/(sum(h)-h(1)), h1(2:end)*(2:length(h1))'/sum(h1(2:end)));
%.1f%% of unique FR1 reads, %.1f%% of all reads).\n', sum(h)-h(1),100*(sum(h)-h(1))/sum(h), 100*h(2:end)*(2:length(h))'/N);
%fprintf('%5d reads appear more than once, but on the same FR1 replicate (%.1f%% of unique FR1 reads, %.1f%% of all reads).\n', sum(h)-h(1),100*(sum(h)-h(1))/sum(h), 100*h(2:end)*(2:length(h))'/N);
subplot(3,1,2);
bar(h/sum(h1));
title('reads in only one replicate, FR1');

ix = sum(counts>0,2) == 1 & sum(counts(:,3:6),2) >= 1;
h1 = sum(counts(ix,:),2);
h1 = hist(h1, 1:max(h1));
fprintf('\t\t%5d reads (%.1f%%) in FR2, %.1f duplicates on average\n', (sum(h1)-h1(1)), 100*(sum(h1)-h1(1))/(sum(h)-h(1)), h1(2:end)*(2:length(h1))'/sum(h1(2:end)));
%fprintf('%d reads appear more than once, but on the same FR2 replicate (%.1f%% of unique FR2 reads, %.1f%% of all reads).\n', sum(h)-h(1),100*(sum(h)-h(1))/sum(h), 100*h(2:end)*(2:length(h))'/N);
subplot(3,1,3);
bar(h1/sum(h1));
title('reads in only one replicate, FR2');

ix = sum(counts>0,2) > 1;
fprintf('\t%6d (%.1f%%) appear in more than one replicate (%.1f%% of all reads).\n', sum(ix),100*sum(ix)/N_, 100*sum(sum(counts(ix,:)))/N );

iy = sum(counts>0,2) > 1 & sum(counts(:,1:2),2) == 0;
fprintf('\t\tOf those %.1f%% from FR2 only', 100*sum(iy)/sum(ix) );

iy = sum(counts>0,2) > 1 & sum(counts(:,3:6),2) == 0;
fprintf(' and %.1f%% from FR1 only', 100*sum(iy)/sum(ix) );
iy = sum(counts>0,2) > 1 & sum(counts(:,3:6),2) > 0 & sum(counts(:,1:2),2) > 0;
fprintf(' (%.1f%% from both).\n', 100*sum(iy)/sum(ix) );

junk = counts(ix,:); 

fprintf('\t\teach replicate has %.1f reads on average, and total copies per read is %.1f on average.\n', mean(junk(junk>0)), mean(sum(junk,2)));


