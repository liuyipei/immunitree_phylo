% data_variants = variant repertoire, an input to get_big_hist
% data_clones = clones repertoire or keep it empty.
% nTrials = number of subsampling trials.
%
% if data_clones is empty, it compares for every individual each tissue 
% repertoire vs.the combined repertoire of the other 3.  
% 
% is data_clones is not empty, it compares for every (individual,tissue) 
% the clones repertoire vs. the variant repertoire.  
function comp = compare_repertoires(data_variants, data_clones, nTrials)
    if ~exist('data_clones', 'var'), data_clones = []; end
    
    if ~isempty(data_clones)
        comp = compare_each_sample(data_variants, data_clones, nTrials);
    else
        one_vs_all = false; % false means all vs all 
        comp = compare_tissue_samples(data_variants, one_vs_all, nTrials);
    end
end

% aggreagates the measures across all trials.  
% all_comp = the results from all trials (an array of the "comp" structure)
% comp = ?
function comp = aggregate_comps(all_comp, comp)    
    if length(all_comp) == 1, comp = all_comp; return; end
    if nargin<2, comp.p = 0; comp.l = 0; comp.kl = 0; comp.l1 = 0; comp.diff = 0; end
    
    comp.p = mean([all_comp.p]) - comp.diff;
    comp.stdp = std([all_comp.p]);
    
    comp.cntl = sum(comp.l<=[all_comp.l]);
    comp.l = mean([all_comp.l]) - comp.l;
    comp.stdl = std([all_comp.l]);
    
    comp.cntkl = sum(comp.kl>=[all_comp.kl]);
    comp.kl = mean([all_comp.kl]) - comp.kl;    
    comp.stdkl = std([all_comp.kl]);
    
    comp.cntl1 = sum(comp.l1>=[all_comp.l1]);
    comp.l1 = mean([all_comp.l1]) - comp.l1;
    comp.stdl1 = std([all_comp.l1]);
    
    comp.diff = mean([all_comp.diff],2) - comp.diff;
    comp.stddiff = std([all_comp.diff]')';
end


% X is the "big" sample, Y is the "small" one
% it computes various of measures: (see below)
function comp = compare_two_samples(X,Y)
        assert(size(X,2) == 1);
        assert(size(Y,2) == 1);
        X_norm = X+0.01;  
        X_norm = X_norm/sum(X_norm);
        Y_norm = Y+0.01;  
        Y_norm = Y_norm/sum(Y_norm);

% for generating Y from the distribution in X:
% comp.p = p-value for generating something as extreme as Y
% comp.l = average log-likelihood for each element in Y
% comp.stddiff(i) = the standard deviation of Y(i)/N if Y was generated
%                   from the null hypothesis
        [comp.p l comp.stddiff] = discrete_test(X_norm, Y);
        comp.l = l / sum(Y);
        
% KL metric difference between empirical and "real" distribution
        comp.kl = KL(X_norm,Y_norm) + KL(Y_norm,X_norm);
        
% difference between empirical distribution and "real" distribution
        comp.diff = Y_norm-X_norm;

% L1 difference between empirical distribution and "real" distribution
        comp.l1 = sum(abs(comp.diff));
end

% compares for each sample the clone and variant repertoire 
function comp = compare_each_sample(data_variants, data_clones, nTrials)
    hist_clones = reshape(get_big_hist(data_clones, false, 0), [], 40);
    hist_variants = reshape(get_big_hist(data_variants, false, 0), [], 40) ;        

    en = max([data_variants.subject_num]);
    for k=1:en
        X = hist_variants(:,k);        
        Y = hist_clones(:,k);        
        comp(k) = compare_two_samples(X,Y);
        if nTrials > 0
            for t=1:nTrials
                X_ = subsample(X,sum(Y), true);
                comp_(k,t) = compare_two_samples(X,X_);
            end
            newcomp(k) = aggregate_comps(comp_(k,:), comp(k));
        end
    end    
    if nTrials > 0, comp=newcomp; end
end

% returns X, a histogram subsampled from the input histogram (also X)
% k = number of entries in the sub-sample
% with_replacement = boolean true/false
function X = subsample(X,k, with_replacement)
    if ~with_replacement
        N = sum(X);
        ix = randperm(N);
        ix = sort(ix(1:k));
        X_ = [0; cumsum(X)+1];
        X = histc(ix,X_)';    
        assert(X(end) == 0);
        X = X(1:end-1);
    else        
        X_norm = X+0.01;  
        X_norm = X_norm/sum(X_norm);
        logX = log(X_norm);
        X = hist(super_lmnrnd(logX', k), 1:length(X))';        
    end
end

% compares one variant sample (of an individual) against all the other 3
% or
% compares all variant samples against each other (one_vs_all = false)
% nTrials = number of sub sampling trials.
function comp = compare_tissue_samples(data_variants, one_vs_all, nTrials)
    if ~exist('nTrials', 'var') || nTrials ==0 || isempty(nTrials)
        sub_sample = 0;
        nTrials = 1;
    else
        sub_sample = 1;
    end
    
%     if one_vs_all, allcomp = zeros(4,10,nTrials);
%     else allcomp = zeros(40,40,nTrials); end
    
    for t=1:nTrials
        if mod(t, 10) == 0, t, end
        hist_variants = get_big_hist(data_variants, false, sub_sample);            
        hist_variants(:,2,1) = 0;
        if one_vs_all        
            for k=1:10                
                for s=1:4
                    Y = hist_variants(:,s, k);
                    X = sum(hist_variants(:,setdiff(1:4, s), k),2);
                    allcomp(s,k,t) = compare_two_samples(X,Y);        
                end
            end   
        else
            %assert(false);
            fprintf('All vs ALL\n');
            hist_variants = reshape(hist_variants, [], 40);
            for k=1:40
                Y = hist_variants(:,k);
                if k==1, allcomp(40,40,nTrials) = compare_two_samples(Y,Y); end   
                for k_=k+1:40                
                    X = hist_variants(:,k_);
                    allcomp(k,k_,t) = compare_two_samples(X,Y);        
                end
            end   
        end
    end
    
    for k=1:size(allcomp,2)            
        for s=1:size(allcomp,1)
            comp(s,k) = aggregate_comps(allcomp(s,k,:));
        end
    end
    
end

function res = binomial(Y)
    res = gammaln(sum(Y)+1) - sum(gammaln(Y+1));
end

function dist = KL(P,Q)
    Q = Q/sum(Q);
    P = P/sum(P);
    temp =  P.*log(P./Q);
    temp(isnan(temp))=0;% resolving the case when P(i)==0
    dist = sum(temp);
end

function test
%% get data,data_clones, data_variants
load ~/scratch/data/lymph/t5.fa_results/combined_clone_repertoire.mat

%%  SUPP FIGURE 1b: L1 divergence between source distributions (with sub sampling)
close all
comp = compare_repertoires(data_variants, [], 1000);
save_figures = 0;
fld = 'l1';
comp(2).(fld) = 0;
if isfield(comp, fld)
    comp(2).(['std' fld]) = 0;
    h = errorb(reshape([comp.l1], 4, 10)', reshape([comp.stdl1], 4, 10)');
else
    h = bar(reshape([comp.l1], 4, 10)');
end
color_the_bars(h);

xlabel('patients');
title('L1 distance between each source VJ normalized repertoire and the rest');
legend({'blood', 'spleen', 'mes', 'med'});
if save_figures
    output_dir = sprintf('~/www/research/donors/nonclones/')
    saveas(gcf, sprintf('%s/%s', output_dir, 'VJ_l1_diff_one_src_vs_all.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'VJ_l1_diff_one_src_vs_all.fig'));
end


%%  Blood vs Rest, no need to do sub-sampling 
labels = {'blood', 'spleen', 'mes', 'med'};
comp = compare_repertoires(data_variants, [], 0);
close all
s = 4;
diff = [comp.diff];
diff = diff(:,s:4:end);
stddiff = [comp.stddiff];
stddiff = stddiff(:,s:4:end);
for k=1:10    
    subplot(10,1,k);
    ix = find(diff(:,k).^2 > max(0.01^2, (2*stddiff(:,k)).^2));
    if 0
        imagesc(reshape(diff(:,k), 6, 57), [-0.05 0.05]);    
        hold on
    else
        % alternative repertoire
        heatmap_to_scatter(reshape(diff(:,k), 6, 57));
    end

    % probability difference more than two std-dev
    [J V] = ind2sub([6 57], ix);
    plot(V, J, 'cx', 'linewidth', 2);
end

subplot(10,1,1);
title(sprintf('diff between normalized repertoires of %s vs the rest combined', labels{s}));

if save_figures
    output_dir = sprintf('~/www/research/donors/nonclones/')
    saveas(gcf, sprintf('%s/%s', output_dir, 'VJ_diff_blood_vs_all.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'VJ_diff_blood_vs_all.fig'));
end


%%  Supp Figure 1D: Clones vs variants.  Comparing clones to simulated clones
close all
comp = compare_repertoires(data_variants, data_clones, 1000);
save_figures = 0;
%
fld = 'l';
comp(2).(fld) = 0;
%[comp(1:4:40).(fld)] = deal(0);
if isfield(comp, ['std' fld])
    comp(2).(['std' fld]) = 0;
%    [comp(1:4:40).(['std' fld])] = deal(0);
    h = errorb(reshape([comp.(fld)], 4, 10)', reshape([comp.(['std' fld])], 4, 10)');
else
    h = bar(reshape([comp.(fld)], 4, 10)');
end
color_the_bars(h);
xlabel('donors');
title(sprintf('E[ %s(simulated clones from variants || variants) - %s(actual clones || variants) ]',...
    fld, fld));
legend({'blood', 'spleen', 'mes', 'med'});
if save_figures
    output_dir = sprintf('~/www/research/donors/combined/')
    filename = 'll_diff_between_mean_variant_repertoire_and_clones_with_error_bars';
    saveas(gcf, sprintf('%s/%s%s', output_dir, filename, '.jpg'));
    saveas(gcf, sprintf('%s/fig/%s%s', output_dir, filename, '.fig'));
end


%%  Clones Vs. Variants, this time focusing on the distribution difference
rep = load_repertoire('ihmmune_collapsed');
labels = {'blood', 'spleen', 'mes', 'med'};
comp = compare_repertoires(data_variants, data_clones, 0);
save_figures = 0;
close all
for s=1:4
    figure
    set(gcf, 'position', [746 76 814 1421]);
    diff = [comp.diff];
    stddiff = [comp.stddiff];
    diff = diff(:,s:4:end);
    stddiff = stddiff(:,s:4:end);
    for k=1:10    
%        subplot(7,1,2*s-1); %k
        subplot(10,1,k);
        if 0
            imagesc(reshape(diff(:,k), 6, 57), [-0.05 0.05]);    
            hold on
        else
            % alternative repertoire
            heatmap_to_scatter(reshape(diff(:,k), 6, 57));%, rep);
        end

        % probability difference more than two std-dev
        ix = find(diff(:,k).^2 > max(0.01^2, (2*stddiff(:,k)).^2));
        [J V] = ind2sub([6 57], ix);
        plot(V, J, 'cx', 'linewidth', 2);

    end
    
    subplot(10,1,1);
    title(sprintf('diff between normalized repertoires of %s clones vs %s variants', ...
        labels{s}, labels{s}));

    if save_figures
        output_dir = sprintf('~/www/research/donors/combined/');
        filename = sprintf('VJ_diff_%s_clones_vs_variants',labels{s});
        saveas(gcf, sprintf('%s/%s.jpg', output_dir, filename));
        saveas(gcf, sprintf('%s/fig/%s.fig', output_dir, filename));
    end

end

%%  Clones Vs. Variants, aggregating all patients
ix = data_clones.subject_num <=40;
data_clones.subject_num(ix) = mod(data_clones.subject_num(ix)-1, 4)+1;
data_variants.subject_num = mod(data_variants.subject_num-1, 4)+1;
comp = compare_repertoires(data_variants, data_clones, 0);
save_figures = 0;
close all
figure
set(gcf, 'position', [746 76 814 1421]);
diff = [comp.diff];
stddiff = [comp.stddiff];
for s=1:4
    subplot(11,1,(3*s-2):(3*s-1));    
    X = reshape(diff(:,s), 6, 57);
    X = [X sum(X,2)/57; sum(X,1)/6 0];
    heatmap_to_scatter(X, rep);
    plot([0 59], [6.5 6.5]);
    plot([57.5 57.5], [0 8]);
    
    % probability difference more than two std-dev
    ix = find(diff(:,s).^2 > max(0.01^2, (2*stddiff(:,s)).^2));
    [J V] = ind2sub([6 57], ix);
    plot(V, J, 'cx', 'linewidth', 2);
    
end
if save_figures
    output_dir = sprintf('~/www/research/donors/combined/');
    filename = sprintf('VJ_diff_aggregated_clones_vs_variants',labels{s});
    saveas(gcf, sprintf('%s/%s.jpg', output_dir, filename));
    saveas(gcf, sprintf('%s/fig/%s.fig', output_dir, filename));
end



%%  Just report the differences per tissue
hist_clones = reshape(get_big_hist(data_clones, false, 0), [], 40);
hist_variants = reshape(get_big_hist(data_variants, false, 0), [], 40) ;        

rep = load_repertoire('ihmmune_collapsed');
fout = 1; %fopen('~/www/research/donors/clones/VJ_frequency_report.txt', 'w');
for s=1:4
    % get repertoire
    cur_sample = s:4:40;

    Y = sum(hist_variants(:,cur_sample),2);
    Y = Y/sum(Y);
    X = sum(hist_clones(:,cur_sample),2);
    X = X/sum(X);
    D = X-Y;
    [~, o] = sort(abs(D), 'descend');
    D = D(o);
    [j v] = ind2sub([6 57], o);
    
    fprintf(fout, '\nsource = %s\n', labels{s});
    fprintf(fout, '========================================\n');
    for i=1:length(v)
        if abs(D(i)) < 0.002, break; end
        fprintf(fout, '%s  \t%s  \t%.2f%% (%.2f vs %.2f)\n', rep.V(v(i)).Header, ...
            rep.J(j(i)).Header, 100*D(i), 100*X(o(i)), 100*Y(o(i)));
    end
end
%fclose(fout);

%%  Supp Figure 1X: Every sample against every sample
close all
comp = compare_repertoires(data_variants, [], 1);
save_figures = 0;
fld = 'l1';
%
Z = zeros(40,40);
for i=1:40,
    for j=1:40
        if ~isempty(comp(i,j).l1)
            Z(i,j) = comp(i,j).l1;
        end
    end
end
Z_ = Z;
%Z = reshape([old_comp.l1], 40, 40);
Z(2,:) = 0;
Z(:,2) = 0;

Z(isnan(Z)) = 0;
Z = Z + Z';

% ix = [2:4:40; 3:4:40; 4:4:40];
% ix = ix(:);
% Z = Z(ix,ix);

output_dir = sprintf('~/www/research/donors/combined/');
figure(1);
imagesc(Z, [min(Z(Z>0)) max(Z(:))]);
plot_class_lines(1:4:size(Z,1)+1, true, 'm', true);
plot_class_lines(1:4:size(Z,2)+1, false, 'm', true);
title('L1 distance between each two of the 39 samples');
xlabel('ordered by donor,tissue');
colorbar 
if save_figures
    output_dir = sprintf('~/www/research/donors/combined/');
    filename = 'L1_donors';
    saveas(gcf, sprintf('%s/%s.jpg', output_dir, filename));
    saveas(gcf, sprintf('%s/fig/%s.fig', output_dir, filename));
end

figure(2)
Z = reshape(Z,4,10,4,10);
Z = permute(Z, [2 1 4 3]);
Z = reshape(Z, 40,40);
imagesc(Z, [min(Z(Z>0)) max(Z(:))]);
plot_class_lines(1:10:size(Z,1)+1, true, 'm', true);
plot_class_lines(1:10:size(Z,2)+1, false, 'm', true);
title('L1 distance between each two of the 39 samples');
xlabel('ordered by tissue, donor');
colorbar 
if save_figures
    output_dir = sprintf('~/www/research/donors/combined/');
    filename = 'L1_tissues';
    saveas(gcf, sprintf('%s/%s.jpg', output_dir, filename));
    saveas(gcf, sprintf('%s/fig/%s.fig', output_dir, filename));
end

end



function obselete
%% Show difference between distributions
Z1 = plot_patient_repertoire(get_big_hist(data_clones, true, 0), {rep.V.Header}, {rep.J.Header});
Z2 = plot_patient_repertoire(get_big_hist(data_variants, true, 0), {rep.V.Header}, {rep.J.Header});
figure;
imagesc(Z1-Z2, [-0.02 0.02]);
colormap default;
plot_class_lines(1:length(rep.V):size(Z1,2)+1, false, 'm', true);
plot_class_lines(1:length(rep.J):size(Z1,1)+1, true, 'm', true);
colorbar
rmse = sqrt(mean((Z1(:)-Z2(:)).^2));
title(sprintf('difference between variant repertoire and clone repertoire\nRMSE = %f', rmse));
pos = get(gcf, 'position');
pos(3:4) = [1024 768];
set(gcf, 'position', pos);
if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'diff_repertoires_variants_clones.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'diff_repertoires_variants_clones.fig'));
end

%% Show in which sample clones distribution different than variants
% likelihood 
[cnt l l_] = compare_repertoires(data_variants, data_clones);

% get how many clones were in each of the 40 samples (except donor 1 spleen)
big_hist = get_big_hist(data_clones, false); 
big_hist = reshape(big_hist, [], 40);
sumClones = sum(big_hist,1)';
%sumClones(2) = 1;

% average clone log-likelihood, per sample
ll = l./sumClones; 

% expected log_likelihood clone, based on variants distribution, per sample
ll_ = l_./sumClones(:,ones(1,size(l_,2)));
ll_mean = mean(ll_,2); %log(mean(exp(ll_),2));

% draw figure
imagesc(reshape(ll-ll_mean, 4, 10)');
title(sprintf('clone empirical log-likelihood minus \nexpected log-likelihood, per sample'));
xlabel('sources');
ylabel('donors');
colorbar('location', 'southoutside')

if save_figures
    saveas(gcf, sprintf('%s/%s', output_dir, 'll_diff_between_mean_variant_repertoire_and_clones.jpg'));
    saveas(gcf, sprintf('%s/fig/%s', output_dir, 'll_diff_between_mean_variant_repertoire_and_clones.fig'));
end




%%
function separate_datas_obselete
% get data_clones
data_ = data;
flds = fields(data_);
for i=1:length(flds)
    fld = flds{i};
    data_.(fld) = data_.(fld)(end-8000+1:end,:);
end
% data_variants
data__ = data;
flds = fields(data__);
for i=1:length(flds)
    fld = flds{i};
    data__.(fld) = data__.(fld)(1:end-8000,:);
end

end



% if ~isempty(data_clones)
%     hist_clones = reshape(get_big_hist(data_clones, false, 0), [], 40);
%     hist_variants = reshape(get_big_hist(data_variants, false, 0), [], 40);
%         
% 
%     nTrials = 1000;  % number of proto-trials for each sample
%     l_ = zeros(40, nTrials);
%     l = zeros(40,1);
%     cnt = l;
%     for k=1:40
%         % variants data, turn each sample into a distribution X over the
%         % V,J combinations
%         X = hist_variants(:,k) + 0.01;
%         X = X(:) / sum(X(:));
% %        X = 0.99*X(:) / sum(X(:))+ 0.01/length(X);
%         
%         % clones data, bin each sample into a histogram Y storing VJ repertoire
%         Y = hist_clones(:,k);        
%         [cnt(k) l(k)] = discrete_test(X, Y);
% 
%         if 0
%            logX = log(X);
%            l(k) = binomial(Y) + sum(Y.*logX); % P(Y | X)
% 
%             for trial = 1:nTrials
%                 % generate sum(Y) samples from X, compute likelihood of
%                 % resulting repertoire
%                 Y_ = hist(super_lmnrnd(logX', sum(Y)), 1:length(X))';
%                 l_(k, trial) = binomial(Y_) + sum(Y_.*logX);
%             end
%         end
%         
%     end
%     
%     if 0
%         suml = sum(l); % total likelihood of clone data given variant data
% 
%         cnt = 0;
%         for k=1:1e5  % real trials
%             % choose for each sample one of its 1000 proto-trials
%             rnd = ceil(nTrials*rand(40,1));
%             t = (1:40)' + 40*(rnd-1);
% 
%             % count how often the joint likelihood of the sampled repertoires
%             % was worse than the joint likelihood of the observed clone
%             % repertoires        
%             cnt = cnt + sum(l_(t))<suml;   
%         end
%     end
%     
% elseif 1 % no data_clones, only variants
%     % compare each tissue the the other 3 tissues combined
%     % output: 
%     %   cnt = p-value for tissue s drawn from the other tissues distribution
%     %   l =  KL distance between the distribution
%     %   l_ = L1 distance
%     cnt = zeros(4,10);
%     l = zeros(4,10);
%     l_ = zeros(4,10);
%     hist_variants = get_big_hist(data_variants, false, 0) + 0.01;
%     for s=1:4                
%         for k=1:10
%             Y = hist_variants(:,s, k);% + 0.01;
%             X = sum(hist_variants(:,[1:s-1 s+1:4], k),2);% + 0.01;
%             X = X/sum(X);
%             [cnt(s,k) l(s,k)] = discrete_test(X, Y);
%             l(s,k) = KL(X,Y) + KL(Y,X);
%             Y = Y/sum(Y);
%             l_(s,k) = sum(abs(Y-X));
%         end
%     end   
% else % no data_clones, only variants
%     % compare each tissue with each tissue
%     cnt = zeros(4,4,10);
%     l = zeros(4,4,10);
%     l_ = zeros(4,4,10);
%     hist_variants = get_big_hist(data_variants, false, 0);
%     for s=1:4                
%         for s_=s+1:4
%         s, s_,
%         for k=1:10
%             Y = hist_variants(:,s, k) + 0.01;
%             X = hist_variants(:,s_, k) + 0.01;        
%             X = X/sum(X);
%             [cnt(s,s_,k) l(s,s_,k)] = discrete_test(X, Y);
%             l(s,s_,k) = KL(X,Y) + KL(Y,X);
%             Y = Y/sum(Y);
%             l_(s,s_,k) = sum(abs(Y-X));
%         end
%     end
%     end
%     
% end
%     
% end
% 
end