%%  Load data and repertoire
clear;
fasta_file_name = '/afs/cs/u/joni/scratch/data/lymph/t5.fa';
load([fasta_file_name '.mat']);   

%% which genes do not even have a single read mapped to them
rep = load_repertoire('ihmmune', true);
[redundant no_show] = find_redundant_Vs(data.ihmmune(:,[1 4:7]), rep);
fprintf('%d V genes are just duplicates of other V genes ', length(redundant));
fprintf('%d V genes do not show up at all\n', length(no_show));
keep = true(1,length(rep.V));
keep(no_show) = false;

%% Do not include pseudo genes  (4624 reads are mapped to them by ihmmune)
for i=1:length(rep.V)
    if ~isempty(find(rep.V(i).Header == 'P',1))
        fprintf('--> '); 
        keep(i) = false;
    end; 
    fprintf('%s\n', rep.V(i).Header); 
end

%% add 'ihmmmune_collapsed' field to the data and create a new repertoire
[rep newrep data] = collapse_repertoire(rep, keep, data);
%% fix IGHV1-c to start in frame and be in par with IMGT aligne reads
assert(isequal(newrep.V(10).Header,'IGHV1-c'));
assert(length(newrep.V(10).Sequence) == 260);
assert(isequal(newrep.V(10).Sequence(1:13), 'GGAAGTCTGGGGC'));
newrep.V(10).Sequence = ['CAGGTCCAGCTGGTGCAGTCTTGGGCTGAGGTGA' newrep.V(10).Sequence];
fprintf('Added 34 bases to IGHV1-c.\n');

%% write new repertoire to disk
rep_dir = '/afs/cs/u/joni/scratch/data/VDJrepertoire';
system(sprintf('mv %s/%s/*.fa %s/%s/bu/', rep_dir, newrep.code));
fastawrite(sprintf('%s/%s/V.fa', rep_dir, newrep.code), newrep.V);
fastawrite(sprintf('%s/%s/D.fa', rep_dir, newrep.code), newrep.D);
fastawrite(sprintf('%s/%s/J.fa', rep_dir, newrep.code), newrep.J);
fprintf('Wrote new repertoire (old one backed up).\n');

%%  unique each read
N = length(data.reads);
iv = [1 4:7];
vs = zeros(N, length(iv));
for i=1:N
    if mod(i,10000) == 0, fprintf('%d ',i); end
    x = unique(data.ihmmune_collapsed(i,iv));
    if x(1) == 0
        vs(i,1:length(x)-1) = x(2:end);
    else
        vs(i,1:length(x)) = x;
    end               
end
fprintf('\n');
data.ihmmune_collapsed(:,iv) = vs;
%%
% save([fasta_file_name '.mat'], 'data');






%% %% %% %% Investigate the difference between the V alleles %% %% %% %%

%% show the distance between the different V alleles

dist = seqpdist({rep.V.Sequence},'Method','jukes-cantor',...
  'Indels','score',...
  'Alphabet', 'NT',...
  'PairwiseAlignment',true);

%% phylo tree of all 330 genes
Tree = linkage(dist, 'complete');
[~,~, order] = dendrogram(Tree, 0);
dist_ = squareform(dist);
imagesc(dist_(order,order), [0 1])

truncated_dist_ = squareform(truncated_dist);
figure;
imagesc(truncated_dist_(order,order), [0 1])

%% phylo tree of all 330 genes tuncated 132 bases
truncated_Vs = cellfun(@(x) x(133:end), {rep.V.Sequence}, 'uniformoutput', false);
truncated_dist = seqpdist(truncated_Vs,'Method','jukes-cantor',...
  'Indels','score',...
  'Alphabet', 'NT',...
  'PairwiseAlignment',true);
truncated_dist_ = squareform(truncated_dist);
truncated_dist_ = squareform(truncated_dist_(keep,keep));
truncated_tree = seqlinkage(truncated_dist, 'complete', {rep.V.Header});
view(truncated_tree);


%% phylo tree of 270 genes (removing those no one maps to)
rep = load_repertoire('ihmmune', true);
rep.V = rep.V(keep);
[~,o] = sort({rep.V.Sequence});
rep.V = rep.V(o);
dist = seqpdist({rep.V.Sequence},'Method','jukes-cantor',...
  'Indels','score',...
  'Alphabet', 'NT',...
  'PairwiseAlignment',true);
dist_ = squareform(dist);

%%%%%%%%%%%%%%%%%%%%%%%%   Done        %%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%% Look at the allele groups of the uncollapsed repertoire %%%%%%%%

%%  Use current uncollapsed repertoire, but order genes by allele-groupes
[rep newrep data] = collapse_repertoire(rep, keep, data);
[sorted_V oV] = sort(rep.map_V);
[sorted_J oJ] = sort(rep.map_J);
oV = oV(sorted_V<max(sorted_V));  % discard Vs with 0 reads mapped to them

%%
data = correct_PCR_artifacts(data); % preprocess data to account for "jackpots"

%%  Plot uncollapsed repertoire over all reads, alleles grouped together
figure;
[~, ~, big_hist] = collect_VJs(data, get_reads(data,1:40), 'ihmmune', 0);%, collapsed)
plot_patient_repertoire(big_hist(oV,oJ), {rep.V(oV).Header}, {rep.J(oJ).Header});
plot_class_lines(sorted_V, false, 'm', false);
plot_class_lines(sorted_J, true, 'm', false);
set(gca, 'position', [0.11 0.11 0.775 0.515]);


%%  plot uncollapsed repertoire over all reads *per patient* (aggregating all sources)
figure;
for k = 1:10
    subplot(10,1,k)
    loc = 1:4;
    cur_sample = (k-1)*4+loc;    
    [~, ~, big_hist] = collect_VJs(data, get_reads(data,cur_sample), 'ihmmune', 0);
    if k == 1    
        plot_patient_repertoire(big_hist(oV,oJ), {rep.V(oV).Header}, {rep.J(oJ).Header});
    else
        plot_patient_repertoire(big_hist(oV,oJ), [], {rep.J(oJ).Header});        
    end
    plot_class_lines(sorted_V, false, 'm', false);
    plot_class_lines(sorted_J, true, 'm', false);
end

%%  plot separate V and J repertoires over all reads per patient
sumV = zeros(10,length(oV)); % sum per patient
sumJ = zeros(length(oJ),10); % sum per patient
for k = 1:10
    loc = 1:4;
    cur_sample = (k-1)*4+loc;    
    [~, ~, big_hist] = collect_VJs(data, get_reads(data,cur_sample), 'ihmmune', 0);
    sumV(k,:) = sum(big_hist(oV,oJ)');
    sumJ(:,k) = sum(big_hist(oV,oJ)',2);
end
 
figure;
plot_patient_repertoire(sumV', {rep.V(oV).Header}, []);
plot_class_lines(sorted_V, false, 'm');
set(gca, 'position', [0.11 0.11 0.775 0.515]);


figure;
plot_patient_repertoire(sumJ', [], {rep.J(oJ).Header});
plot_class_lines(sorted_J, true, 'm', false);
set(gca, 'XTick', 1:10)
title('J genes for each patient');

%%  Show the length of each V gene
    figure;
    subplot(4,1,4);    
    bar(sum(sumV));
    plot_class_lines(sorted_V, false, 'm');
%    plot_patient_repertoire(sumV', {rep.V(oV).Header}, []);
    xlim([0 270]);

    %figure
    subplot(4,1,1:3);
    lens = cellfun(@length, {rep.V(oV).Sequence});      
    bar(lens);
    plot_class_lines(sorted_V, false, 'm');
    ylim([150 350]);
    xlim([0 270]);
    h = text(1:270, 300+50*mod(1:270,2), ...
        cellfun(@(x) x(4:end), {rep.V(oV).Header}, 'uniformoutput', false));
    set(h, 'rotation', 90);
    set(h, 'fontsize', 8);
    %title('length of V genes in repertoire')    

%%%%%%%%%%%%%%%%%%%%%%%%   Done        %%%%%%%%%%%%%%%%%%%%%%%%%
    

%% remember original uncollapsed repertoire
ihmmune_bu = data.ihmmune;  % back up uncollapsed alignment

%%  collapse repertoire! (change data.ihmmune)
%rep = collapse_repertoire(rep, keep);
ihmmune = ihmmune_bu(:,[1 4:7]);
ix = ihmmune > 0;
ihmmune(ix) = rep.map_V(ihmmune(ix));
data.ihmmune(:, [1 4:7]) = ihmmune;
ix = ihmmune_bu(:,3) > 0;
data.ihmmune(ix,3) = rep.map_J(ihmmune_bu(ix,3));
