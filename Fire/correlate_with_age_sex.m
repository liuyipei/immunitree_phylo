lbl = [];
lbl.num = 2;
lbl.names = {'0','1'};
lbl.str = 'class';
lbl.patients = [repmat([1 2 3 4 5], 100,1) repmat([1 2 3 4 5], 100,1)];
lbl.patients = lbl.patients(:);
lbl.l = [ones(500,1); 2*ones(500,1)];
y = 100 + randn(1000,1);



%%
fasta_file_name = '/afs/cs/u/joni/scratch/data/lymph/t5.fa';
load([fasta_file_name '.mat']);

%% playing with random effect and correlating to subject demographies

X = clone_intra_read(:,5);
X = dist_to_germ; %';
x = X;
L = 10;
for s=-1:5    
    if s == -1
        ix = (labels.l >1); %== s);
    elseif s == 0
        ix = true(length(x),1); 
    else
        ix = (labels.l == s);
    end
    
    T = zeros(1,L); v = zeros(1,L);
    for k=1:10   
        iy = labels.patients==k & ix;
        T(k) = mean(x(iy));
        v(k) = var(x(iy));
    end
    fprintf('\ns=%d (total %d):\n',s, sum(ix));
    fprintf('\t%.2f', T);
    fprintf('\n');
end

%%
[p,z,T_star,v_star] = random_effect(T, v, true);
c = corr(T', data.patients(:,2));
[~, p1] = ttest2(T(data.patients(:,2)<62), T(data.patients(:,2)==62));
[~, p2] = ttest2(T(data.patients(:,1)==0), T(data.patients(:,1)==1));
fprintf('corr = %.2f, pval1 = %0.2f, pval2 = %0.2f\n', c,p1,p2);

%end
%% playing with ANOVA

X = clone_intra_read(:,5);
ix = 1:N;% (labels.patients~=1);
[p,t,stats,terms] = anovan(X(ix),{labels.l(ix) labels.patients(ix)},'varnames', {'source', 'patient'}, 'model', 'linear', 'random', 2);
%multcompare(stats)



%% ANOVA on repertoire with/without considering source
experiment_id = 9;

big_hist_ = reshape(big_hist, length(rep.J), length(rep.V), 4, 10);
% big_hist_ = big_hist_(:,:,2:4,:);
% big_hist_ = big_hist_(:,:,:,2:10);
[J V S P] = ind2sub(size(big_hist_), 1:numel(big_hist_));
factors = [0 0 0 0 ; 1 0 1 0; 1 0 0 1; 0 1 1 0; 0 1 0 1; 1 1 0 0];
    
[p,t,stats,terms] = anovan(big_hist_(:),{J V S P},'varnames', {'J', 'V', 'source', 'patient'}, 'model', factors);%,'random', 4);
fprintf('rmse = %.3g\n', sqrt(mean(stats.resid.^2)));
% for dim=1:4
%     subplot(4,1,dim);
%     multcompare(stats, 'dimension', dim);
%     saveas(gcf, sprintf('ANOVA_repertoire_experiment_%d_dim_%d.jpg', experiment_id, dim));
% end

% Analysis of Variance
%   Source	Sum Sq.     d.f.	Mean Sq.	F	Prob>F
%   V       1617901     56      28891       26.7376	0
%   J       170171  	5       34034       8.0459	0.00083113
%   source	31262       2       15631       2.0724	0.17258
%   patient	172813      9       19201       14.9791	0
%   Error	11007482	10187	1080		
%   Total	12999630	10259			
% Constrained (Type III) sums of squares.


%%  Reconstruction based on the ANOVA decomposition
big_hist__  = zeros(max(J), max(V), max(S), max(P));
all = {1:max(J), 1:max(V), 1:max(S), 1:max(P)};
for i=1:size(stats.vars,1)
    ix = all;
    for j=1:4
        if stats.vars(i,j) ~= 0
            ix{j} = stats.vars(i,j);
        end
    end       
    big_hist__(ix{1},ix{2},ix{3},ix{4}) = big_hist__(ix{1},ix{2},ix{3},ix{4}) + stats.coeffs(i);    
end

figure;
Z_ = plot_patient_repertoire(reshape(big_hist__, [], 4, 10), rep.V, rep.J, 'reconstruction');
fprintf('rmse = %.3g\n', sqrt(mean((Z(:)-Z_(:)).^2)));




%%  VJ not independent ANOVA
big_hist_ = reshape(big_hist, [], 4, 10);
big_hist_ = big_hist_(:,2:4,:);
[VJ S P] = ind2sub(size(big_hist_), 1:numel(big_hist_));
[p,t,stats,terms] = anovan(big_hist_(:),{VJ S P},'varnames', {'VJ', 'source', 'patient'}, 'model', 'linear', 'random', 3);

% Analysis of Variance
%   Source	Sum Sq.	d.f.	Mean Sq.	F	Prob>F
%   VJ	9531273.195	341	27950.9478	84.8303	0
%   source	31262.1976	2	15631.0988	2.7054	0.11612
%   patient	172813.0755	9	19201.4528	58.0437	0
%   Error	3264281.6418	9907	329.4924		
%   Total	12999630.1098	10259			
% Constrained (Type III) sums of squares.


%% chi2 to tell if V and J are independent for each repertoire
p = zeros(10,4);
big_hist = reshape(big_hist, 6,57, []);
for i=[1 3:40]
    X = big_hist(:,:,i);
    p(i) = chi2independence(X);
end


%% For each V, see if there is a source where it is more/less likely
rep = load_repertoire('ihmmune_collapsed');
close all
big_hist_norm = get_big_hist(data, true);
big_hist_norm = reshape(big_hist_norm, 6, 57, 4, 10);
X = sum(big_hist_norm, 1);
VS = zeros(57,4);
for v=1:57
    for s=1:4
        x = X(1,v, s,:);
        y = X(1,v, [1:s-1 s+1:4],:);
        x = x(:); y = y(:);
        % ttest x and y
        [~, VS(v,s)] = ttest2(x, y);
        if VS(v,s)<0.05/(57*4)
            v,s,x,y
            figure;
            hx = hist(x, 0:0.001:0.1);
            hy = hist(y, 0:0.001:0.1);
            bar(0:0.001:0.1, hx, 'r');
            hold on
            bar(0:0.001:0.1, hy, 'g');
        end
    end    
end

%% For each J, see if there is a source where it is more/less likely
close all
big_hist_norm = get_big_hist(data, true);
big_hist_norm = reshape(big_hist_norm, 6, 57, 4, 10);
X = sum(big_hist_norm, 2);
JS = zeros(6,4);
for j=1:6
    for s=1:4
        x = X(j,1, s,:);
        y = X(j,1, [1:s-1 s+1:4],:);
        x = x(:); y = y(:);
        % ttest x and y
        [~, JS(j,s)] = ttest2(x, y);
        if JS(j,s)<0.05%/(6*4)
            j,s,x,y
            figure;
            hx = hist(x, 0:0.01:0.9);
            hy = hist(y, 0:0.01:0.9);
            bar(0:0.01:0.9, hy, 'g');
            hold on
            bar(0:0.01:0.9, hx, 'r');
        end
    end    
end

