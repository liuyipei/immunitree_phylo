function analyze_repertoire

clear;
cd ~/JVL/src/phylo/Fire
addpath('.');
addpath('../../phylo/');
addpath('../../DPtrees/');
addpath('../../VDJ/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));
output_dir = sprintf('~/www/research/donors/combined');

%%  Load data and repertoire
fasta_file_name = '/afs/cs/u/joni/scratch/data/lymph/t5.fa';
load([fasta_file_name '.mat']);      
data = correct_PCR_artifacts(data); % preprocess data to account for "jackpots"
rep = load_repertoire('ihmmune_collapsed');

%%  Plot collapsed repertoire over all reads
figure;
[~, ~, big_hist] = collect_VJs(data, get_reads(data, 1:40), 'ihmmune_collapsed', 0)
plot_patient_repertoire(big_hist, {rep.V.Header}, {rep.J.Header});
set(gca, 'position', [0.11 0.11 0.775 0.515]);
total_entries = sum(big_hist(:));
total_reads = sum(   (data.ihmmune(:,1) ~= 0) &  (data.ihmmune(:,3) ~= 0));
fprintf('total entries in histogram: %d\n', total_entries);
fprintf('total reads with VJ alignment: %d\n', total_reads);
fprintf('Hence %d double counts.\n', total_entries - total_reads);

%%  plot collapsed repertoire over all reads *per patient* (sources aggregated)
figure;
for k = 1:10
    subplot(10,1,k)
    loc = 1:4;
    cur_sample = (k-1)*4+loc;    
    [~, ~, big_hist] = collect_VJs(data, get_reads(data,cur_sample), 'ihmmune_collapsed', 0);%, collapsed)
    if k == 1
        plot_patient_repertoire(big_hist, {rep.V.Header}, {rep.J.Header}, sprintf('patient %d', k));
    else
        plot_patient_repertoire(big_hist, [], {rep.J.Header}, sprintf('patient %d', k));
    end
end
pos = get(gcf, 'position');
pos(3:4) = [1024 768];
set(gcf, 'position', pos);
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'patient_repertoires.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'patient_repertoires.fig'));
end

%%  Show repertoires per source
figure;
source = {'blood', 'spleen', 'LN1', 'LN2'};
rep = load_repertoire('ihmmune_collapsed');
for s=[1 2 3 4]
    subplot(4,1,s);
    [~, ~, big_hist] = collect_VJs(data, get_reads(data,s:4:40), 'ihmmune_collapsed', 0);
    if (s == 1)
        plot_patient_repertoire(big_hist, {rep.V.Header}, {rep.J.Header}, source{s});
    else
        plot_patient_repertoire(big_hist, [], {rep.J.Header}, source{s});
    end

end


%%  Full repertouire per patient and per source
rep = load_repertoire('ihmmune_collapsed');
big_hist = get_big_hist(data, true, 0);
figure;
Z = plot_patient_repertoire(big_hist, {rep.V.Header}, {rep.J.Header}, ...
    'collapsed repertoire for all patients/sources, normalized per sample');
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'all_40_samples.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'all_40_samples.fig'));
end

%% %%%%%%%%%     Factor analysis on matrix    %%%%%%%%%%%%%%%%%%%%%%%%%%%
deep_factorization = true;
mean_factorization = 1;
[rmse a b big_hist_ delta] = factorize_big_hist(big_hist, rep, deep_factorization, mean_factorization );
figure;
Z_ = plot_patient_repertoire(big_hist_, rep.V, rep.J, ...
    'reconstructed repertoire for all patients/sources, normalized per sample');
rmse*1000
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'factor_analysis_reconstruction.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'factor_analysis_reconstruction.fig'));
end
%% Show residue between repertoire constructed from factors and truth
figure;
imagesc(Z-Z_, [-0.02 0.02]);
colormap default;
plot_class_lines(1:length(rep.V):size(Z,2)+1, false, 'm', true);
plot_class_lines(1:length(rep.J):size(Z,1)+1, true, 'm', true);
colorbar
rmse = sqrt(mean((Z(:)-Z_(:)).^2));
title(sprintf('difference between true repertoire and reconstructed repertoire\nRMSE = %f', rmse));
pos = get(gcf, 'position');
pos(3:4) = [1024 768];
set(gcf, 'position', pos);
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'factor_analysis_diff.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'factor_analysis_diff.fig'));
end

%% Another repertoire plot of all samples, but this time with the factors
figure; 
mean_big_hist = mean(reshape(big_hist, [], 40),2);
mean_big_hist = reshape(mean_big_hist(:, ones(1,40)), [], 4, 10);

subplot(5,5, [1:4 6:9 11:14 16:19]);
Z = plot_patient_repertoire(big_hist - mean_big_hist, rep.V, rep.J, ...
    'difference between the sample repertoire and the average repertoire');
colormap default

subplot(5,5,21:24);
imagesc(reshape(a', 6, []));%, [-4 4]);
plot_class_lines(1:length(rep.V):size(Z,2)+1, false, 'm', true);
colorbar('EastOutside');
axis off

subplot(5,5, 5:5:20);
imagesc(reshape(permute(reshape(b', 6, [], 10), [1 3 2]),60,[]));%, [-0.06 0.06]);
plot_class_lines(1:length(rep.J):size(Z,1)+1, true, 'm', true);
colorbar('SouthOutside');
axis off

%% manually adjust location of colorbars
pos = get(gcf, 'position');
pos(3:4) = [1200 800];
set(gcf, 'position', pos);
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'factor_analysis_factors.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'factor_analysis_factors.fig'));
    fprintf('Figures saved.\n');
end



%% Permutation Tests
deep_factorization = true;
mean_factorization = 1;

rep = load_repertoire('ihmmune_collapsed');
nIter = 1000;
permutation_codes = 1:5;
rmses = zeros(nIter,5);
for it = 1:nIter
    it
    big_hist = get_big_hist(data, true, data.dist2germ<=1);   
    for j=permutation_codes
        big_hist_permuted = permute_big_hist(big_hist, j);
        rmses(it,j) = factorize_big_hist(big_hist_permuted, rep, deep_factorization, mean_factorization);        
    end
end

close all
figure;
set(gcf, 'position', [742 1232 818 300]);
for j=2:length(permutation_codes);
    k = permutation_codes(j);
    subplot(length(permutation_codes)-1, 1, k-1);
    y = rmses(:,k)-rmses(:,1);
    maxy = max(abs(y));
    hist(y, -maxy:1e-6:maxy);
    colormap summer
    title(sprintf('%d permutation test %d\np-value = %.3f', nIter, k, sum(y<0)/length(y)));
    xlabel('rmse diff');
    ylabel('number of permutations');
    if save_figures
        saveas(gcf, sprintf('%s/%s', output_dir, 'factor_analysis_permutation1.jpg'));
        saveas(gcf, sprintf('%s/fig/%s', output_dir, 'factor_analysis_permutation1.fig'));
    end
end


%%  Report the best hits for every sample
% for each sample
rep = load_repertoire('ihmmune_collapsed');
fout = 1; %fopen('~/www/research/donors/clones/VJ_frequency_report.txt', 'w');
for d=1:10
    for s=1:4
        % get repertoire
        cur_sample = 4*(d-1)+s;
        [v j X] = collect_VJs(data, get_reads(data, cur_sample), 'ihmmune_collapsed', -20);
        fprintf(fout, '\ndonor = %d source = %d;  %d clones\n', d, s, sum(X(:)));
        fprintf(fout, '========================================\n');
        for i=1:length(v)
            if X(v(i), j(i)) <3, continue; end
            fprintf(fout, '%s  \t%s  \t%.1f%%\n', rep.V(v(i)).Header, rep.J(j(i)).Header, 100*X(v(i), j(i))/sum(X(:)));
        end
    end
end
fclose(fout);



%%  Permutation Test1 - switch sources
% nIter = 150;
% rmse0 = zeros(nIter,1);
% rmse1 = zeros(nIter,1);
% for it = 1:nIter
%     it
%     big_hist = get_big_hist(data, true, 1);   
%     big_hist_permuted = big_hist;
%     % take big_hist and permute the sources for each patient
%     for k=1:10  
%         ix = randperm(4);
%         big_hist_permuted(:,:,k) = big_hist_permuted(:,ix,k);
%     end
%     
%     rmse0(it) = factorize_big_hist(big_hist, rep, deep_factorization,mean_factorization);
%     rmse1(it) = factorize_big_hist(big_hist_permuted, rep, deep_factorization,mean_factorization);
% end
% %%
% close all
% figure;
% y = rmse1-rmse0;
% hist(y, min(y):1e-6:max(y));
% %hist(rmse1, min(rmse1):0.000001:max(rmse1))
% hold on
% lim_y = ylim;
% %hist(rmse0, min(rmse0):0.000001:max(rmse0),'r')
% plot(rmse, lim_y(2)/2, 'xr', 'markersize', 20, 'linewidth', 5);
% colormap summer
% title(sprintf('%d permutation tests over the samples of each individual\np-value = %.3f', nIter, sum(rmse>rmse1)/length(rmse1)));
% xlabel('rmse');
% ylabel('number of permutations');
% set(gcf, 'position', [742 1232 818 300]);
% if save_figures
%     saveas(gcf, sprintf('%s/%s', output_dir, 'factor_analysis_permutation1.jpg'));
%     saveas(gcf, sprintf('%s/fig/%s', output_dir, 'factor_analysis_permutation1.fig'));
% end
% 
% 
% %%  Permutation Test2 - switch non-blood sources
% big_hist_permuted = big_hist;
% nIter = 150;
% rmse2 = zeros(nIter,1);
% for it = 1:nIter
%     it
%     % take big_hist and permute the sources for each patient
%     for k=1:10  
%         ix = 1+randperm(3);
%         big_hist_permuted(:,2:4,k) = big_hist_permuted(:,ix,k);
%     end
%     rmse2(it) = factorize_big_hist(big_hist_permuted, rep, deep_factorization, mean_factorization);
% end
% %%
% figure;
% hist(rmse2, min(rmse2):0.000001:max(rmse2))
% hold on
% lim_y = ylim;
% plot(rmse, lim_y(2)/2, 'xr', 'markersize', 20, 'linewidth', 5);
% colormap summer
% title(sprintf('%d permutation tests over the non-blood samples of each individual\np-value = %.3f', nIter, sum(rmse>rmse2)/length(rmse2)));
% xlabel('rmse');
% ylabel('number of permutations');
% set(gcf, 'position', [742 1232 818 300]);
% if save_figures
%     saveas(gcf, sprintf('%s/%s', output_dir, 'factor_analysis_permutation2.jpg'));
%     saveas(gcf, sprintf('%s/fig/%s', output_dir, 'factor_analysis_permutation2.fig'));
% end
% 
% %%  Permutation Test3 - switch patients
% big_hist_permuted = big_hist;
% nIter = 1000;
% rmse3 = zeros(nIter,1);
% for it = 1:nIter
%     if mod(it,10) == 0, fprintf('.'); end
%     if mod(it,100) == 0, fprintf('\n'); end
%     
%     for k=1:4  
%         ix = randperm(10);
%         big_hist_permuted(:,k,:) = big_hist_permuted(:,k,ix);
%     end
%     rmse3(it) = factorize_big_hist(big_hist_permuted, rep, deep_factorization, mean_factorization);
% end
% %%
% figure;
% hist(rmse3, min(rmse3):0.00001:max(rmse3))
% hold on
% lim_y = ylim;
% plot(rmse, lim_y(2)/2, 'xr', 'markersize', 20, 'linewidth', 5);
% colormap summer
% title(sprintf('%d permutation tests over the patients for each source tissue\np-value = %.3f', nIter, sum(rmse>rmse3)/length(rmse3)));
% xlabel('rmse');
% ylabel('number of permutations');
% set(gcf, 'position', [742 1232 818 300]);
% if save_figures
%     saveas(gcf, sprintf('%s/%s', output_dir, 'factor_analysis_permutation3.jpg'));
%     saveas(gcf, sprintf('%s/fig/%s', output_dir, 'factor_analysis_permutation3.fig'));
% end
% %%  Permutation Test4 - switch VJ combinations
% big_hist_permuted = big_hist;
% nIter = 1000;
% rmse4 = zeros(nIter,1);
% for it = 1:nIter
%     if mod(it,10) == 0, fprintf('.'); end
%     if mod(it,100) == 0, fprintf('\n'); end
%     
%     % take big_hist and permute the sources for each patient
%     ix = randperm(size(big_hist,1));
%     big_hist_permuted = big_hist_permuted(ix,:,:);
%     rmse4(it) = factorize_big_hist(big_hist_permuted, rep, deep_factorization);
% end
% %%
% figure;
% hist(rmse4, min(rmse4):0.00001:max(rmse4))
% hold on
% lim_y = ylim;
% plot(rmse, lim_y(2)/2, 'xr', 'markersize', 20, 'linewidth', 5);
% colormap summer
% title(sprintf('%d permutation tests over the order of VJ combinations\np-value = %.3f', nIter, sum(rmse>rmse4)/length(rmse4)));
% xlabel('rmse');
% ylabel('number of permutations');
% set(gcf, 'position', [742 1232 818 300]);
% if save_figures
%     saveas(gcf, sprintf('%s/%s', output_dir, 'factor_analysis_permutation4.jpg'));
%     saveas(gcf, sprintf('%s/fig/%s', output_dir, 'factor_analysis_permutation4.fig'));
% end
% 
% 








%%   Analysis of out-of-frame reads repertoire?
big_hist = zeros(length(rep.V), length(rep.J), 40);
for k=1:40
    tot = length(get_reads(data,k));
    ix = get_reads(data,k, true, [], true);
    fprintf('%d, %d/%d reads. ', k, length(ix), tot);
    [~, ~, X] = collect_VJs(data, ix, 'ihmmune_collapsed', 0);
    big_hist(:,:,k) = X ./ sum(X(:));
end

big_hist = permute(big_hist, [2 1 3]);
big_hist = reshape(big_hist, [], 4, 10);

% fix patient1/spleen to just be the average of the rest
big_hist(:,2, 1) = mean(big_hist(:,[1 3 4],1), 2);

figure;
Z = plot_patient_repertoire(big_hist, rep.V, rep.J, ...
    'out-of-frame reads repertoires, normalized per sample');


%% %%%%%%%%%     Factor analysis on matrix    %%%%%%%%%%%%%%%%%%%%%%%%%%%
deep_factorization = true;
mean_factorization = 1;
[rmse a b big_hist_ delta] = factorize_big_hist(big_hist, rep, deep_factorization, mean_factorization );
rmse*1000
%%
deep_factorization = false;
mean_factorization = 1;
[rmse a b big_hist_ delta] = factorize_big_hist(big_hist, rep, deep_factorization, mean_factorization );
rmse*1000

big_hist_baseline = mean(reshape(big_hist, size(big_hist,1), []), 2);
big_hist_baseline = big_hist_baseline(:, ones(1,4), ones(1,10));
big_hist_ = big_hist-big_hist_baseline;
rmse = 1000*sqrt(mean( (big_hist(:)-big_hist_baseline(:)).^2))

%%
data_ = data;
ix = find(data.dist2germ <= 1);
flds = fields(data);
for f = 1:length(flds)
    data_.(flds{f}) = data_.(flds{f})(ix,:);
end



end