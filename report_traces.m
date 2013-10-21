function report_traces(a,stats2,stats, canonical_ll)
global codon2aa

figure(2);
subplot(3, 1, 3);
imagesc(a.mut_model.rate_class);
xlim([0 length(a.mut_model.rate_class)]);
title('rate class');
colorbar('East');
nClasses = size(a.mut_model.NT,3);
figure(3); 
row = 7;
for k=1:nClasses
    subplot(nClasses,row,row*(k-1)+(4:5));
    imagesc(a.mut_model.NT(:,:,k)./repmat(sum(a.mut_model.NT(:,:,k)), size(a.mut_model.NT,1),1), [0 0.05]);
    set(gca, 'XTick', 1:4);
    set(gca, 'XTickLabel', ('ACGT')');
    set(gca, 'YTick', 1:4);
    set(gca, 'YTickLabel', ('ACGT')');

    if isfield(a.mut_model, 'AA')
        subplot(nClasses,row,row*(k-1)+(1:3));
        imagesc(a.mut_model.AA(:,:,k)./repmat(sum(a.mut_model.AA(:,:,k)), size(a.mut_model.AA,1),1), [0 0.1]);
        set(gca, 'XTick',1:20);
        set(gca, 'XTickLabel', int2aa(1:20)');
        set(gca, 'YTick',1:20);
        set(gca, 'YTickLabel', int2aa(1:20)');
        colorbar('WestOutside');


        subplot(nClasses,row,row*(k-1)+(6:7));
        Q = make_codon_matrix(a.mut_model.AA(:,:,k), a.mut_model.NT(:,:,k));
        if 1
            imagesc(Q, [0 0.01]);    
        else
            if 0    
                [labels, ord] = sort(codon2aa);
                [labels, I] = unique(labels, 'first');
                imagesc(Q(ord,ord), [0 0.01]);
            else    
                Q_ = zeros(21,21);
                for l=1:21 % child
                    for m=1:21 % parent
                        lx = (codon2aa == l);
                        mx = (codon2aa == m);
                        Q_(l,m) = mean(sum(Q(lx,mx)));
                    end
                end
                imagesc(Q_, [0 0.01]);
                I = (1:21)';  
                labels = I;
            end
            set(gca, 'XTick', I);
            set(gca, 'XTickLabel', int2aa(labels));
            set(gca, 'YTick', I);
            set(gca, 'YTickLabel', int2aa(labels));
        end
    end
end



%%
figure(5);
clf
lim = find(stats2(1,:) == 0, 1);
if isempty(lim), lim = size(stats2,2);
else lim = lim-1;
end

subplot(4,1,1);
plot(stats2(1,1:lim));
title(sprintf('Total %d iterations, with %.2f sec per iteration', lim, max(stats2(end,1:lim))/lim));
stats2(1,min(lim,20):lim);
ylim([min(stats2(1,min(lim,20):lim)), max(stats2(1,min(lim,20):lim))]);
ylabel('log likelihood');


if exist('canonical_ll', 'var')
    subplot(4,1,2);
    plot(canonical_ll);
    title(sprintf('LL of canonical trees'));    
    ylabel('log likelihood');
end

subplot(4,1,3);
junk_ = squeeze(stats(:,1,1:lim));
junk = junk_;
plot((1:numel(junk))/4, junk(:));    
hold on

junk = squeeze(stats(:,2,1:lim));
plot((1:numel(junk))/4, junk(:), 'r');    
plot((1:numel(junk))/4, junk_(:)-junk(:), 'g');
mean_births = mean(junk_(ceil(end/2):end));
mean_deaths = mean(junk(ceil(end/2):end));

junk = squeeze(stats(:,3,1:lim));
plot((1:numel(junk))/4, junk(:), 'm');
title(sprintf('on average: %.1f births and %.1f deaths (%.1f live cells); snapshot time = %.1f', mean_births, mean_deaths, mean_births-mean_deaths, a.F));

legend({'births', 'deaths', 'alive', 'lifetime'}, 'Location', 'NorthWest');
xlabel('iteration');  ylabel('no. nodes');
%ylim([0 max(junk(3,min(lim,40):lim))]);

subplot(4,1,4);
plot(stats2(2,1:lim)./stats2(4,1:lim));
hold on
plot(stats2(3,1:lim)/length(a.t), 'r');
plot(stats2(9,1:lim), 'g');
title('mismatches counts')
xlabel('iterations');
ylim([0 1/length(a.t)*max(1,max(stats2(3,min(lim,10):lim)))]);
legend({'avg mutations per birth', 'avg errors per read', 'RAND index'}, 'Location','SouthWest');

% junk = squeeze(stats(:,3,1:lim));
% plot((1:numel(junk))/4, junk(:)./junk_(:), 'm');


% birth/death rates, similar to births/deaths
% plot(stats2(5,1:lim), 'g');
% plot(stats2(6,1:lim), 'm');

end
