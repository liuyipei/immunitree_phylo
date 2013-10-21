load ~/scratch/data/lymph/t5.fa_results_official/combined_clone_repertoire.mat


%%  CRF Learning
big_hist = get_big_hist(data, false); % unnormalized

scopes{1} = {[1 2 3], [1 2 4], [3 4]};
scopes{2} = {[1 2], [1 4],[1 3], [2 4], [2 3], [3 4]};
scopes{3} = {[1 2 3], [3 4]};
scopes{4} = {[1 2 4], [3 4]};
scopes{5} = {[1 3], [2 3], [3 4]};
scopes{6} = {[1 4], [2 4], [3 4]};
scopes{7} = {[1 2], [3 4]};
scopes{8} = {[1], [2], [3 4]};

big_hist = reshape(big_hist, 6,57, 4,10);
big_hist_ = reshape(big_hist, [], 40);
sum_hist = sum(big_hist_);
big_hist_ = big_hist_ ./ sum_hist(ones(size(big_hist_,1),1), :);
big_hist_norm = reshape(big_hist_, 6,57, 4, 10);

ll_ = [];
phases = [10000 0.01; 100 0.001; 100 0.0001; 1000 0.00001];
rmse1 = zeros(size(phases,1),length(scopes));
rmse2 = zeros(size(phases,1),length(scopes));
for i=1%:length(scopes)
    factors = [];
    for j = 1%size(phases,1)
        [factors reconstruction converged ll] = CRF_learn(big_hist, scopes{i}, phases(j,1),  phases(j,2), factors);
        norm_reconstruction = CRF_expected_counts(size(reconstruction), scopes{i}, factors(1:end-1));
        rmse1(j,i) = sqrt(mean((reconstruction(:)-big_hist(:)).^2))
        rmse2(j,i) = sqrt(mean((norm_reconstruction(:)-big_hist_norm(:)).^2))
        ll_= [ll_; ll];
    end
end


%% Test CRF Learning (delete when done)
scopes = {[1 2 3], [1 2 4], [3 4]};
scopes = {[1 2], [1 4],[1 3], [2 4], [2 3], [3 4]};
big_hist = get_big_hist(data, false); 
big_hist = reshape(big_hist, 6,57, 4,10);
marginals = CRF_marginal_counts(big_hist, scopes);
factors = cell(length(scopes),1);
for i=1:length(factors)-1
    factors{i} = log(marginals{i}+1); % can't allow zeros!
end
factors{end} = log(marginals{end});
        
%%
[factors reconstruction converged ll] = CRF_learn(big_hist, scopes, 20,  0.001, factors_0);


%%  Use MATLAB unconstrained minimization method with the given gradient
f = @(x) delme(x, big_hist, [marginals{1}(:); marginals{2}(:)], repmat(marginals{3}(:)',6*57,1));
options = optimset('fminunc');
options = optimset(options, 'Display','iter-detailed');
options = optimset(options, 'GradObj','on');
options = optimset(options,'LargeScale','off');
[x,fval] = fminunc(f,[factors{1}(:); factors{2}(:)] , options);
%%
norm_reconstruction = CRF_expected_counts(size(big_hist), scopes, factors_0(1:end-1));

%%
%%  Use MATLAB unconstrained minimization method with the given gradient
f = @(x) delme(x, big_hist, ... 
    [marginals{1}(:); ....
    marginals{2}(:); ...
    marginals{3}(:); ...
    marginals{4}(:); ...
    marginals{5}(:)], ...
    repmat(marginals{6}(:)',6*57,1));
options = optimset('fminunc');
options = optimset(options, 'Display','iter-detailed');
options = optimset(options, 'GradObj','on');
options = optimset(options,'LargeScale','off');
[x,fval] = fminunc(f,...
    [factors{1}(:); ...
    factors{2}(:);  ...
    factors{3}(:);  ...
    factors{4}(:);  ...
    factors{5}(:)],  ...
    options);
%%
[J V S D] = size(big_hist);
sz = [J V S D];    
factors_0 = factors;
ix  = [ 0         342         402         426         996        1224];
for i=1:length(factors)-1
    factors_0{i} = reshape(x(ix(i)+1:ix(i+1)), sz(scopes{i}(1)), sz(scopes{i}(2)));
end
%%
[junk g] = f([factors{1}(:); ...
    factors{2}(:);  ...
    factors{3}(:);  ...
    factors{4}(:);  ...
    factors{5}(:)]) 




