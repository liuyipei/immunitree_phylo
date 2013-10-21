% inputs: 
%    param = sicrete distribution (sums to 1)
%    data = count histogram
% output:
%    p = p-value of null hypothesis that data this extreme was generated from param
%    l = log_likelihood of data = log P(data | param)
%    std_dev = an array the size of param.  std_dev(i) is the standard 
%              deviation from param(i) of what we expect data(i)/N to be.
function [p l std_dev p_ l_] = discrete_test(param, data)

assert(numel(param)==numel(data));
assert(abs(sum(param)-1) < 1e-3);

ix = param==0;
data(ix) = [];
param(ix) = [];

N = sum(data);
expected = N*param;

% option 1: GOF
    bins = (1:numel(param))';
    [h,p,stats] = chi2gof(bins, 'frequency', data(:),  'ctrs', bins, 'expected', ...
        expected(:), 'nparam', 0,  'emin', 0);
    std_dev = sqrt(param.*(1-param)/N);

    logParam = log(param);    
    l = binomial(data) + sum(data.*logParam); % log P(data | param)    
    
% option 2: log-likelihood
if 0
    nTrials = 1000;
    dump = zeros(nTrials,length(param));
    l_ = zeros(nTrials,1);
    for i=1:nTrials
        synthetic = hist(super_lmnrnd(logParam', N), 1:length(param))';
        dump(i,:) = synthetic;
        l_(i) = binomial(synthetic) + sum(synthetic.*logParam); % logP(synthetic | param)        
    end
    p_ = sum(l > l_)/nTrials;
    
    std_dev_ = std(dump)/N;
end

end

function res = binomial(Y)
    res = gammaln(sum(Y)+1) - sum(gammaln(Y+1));
end

function test()
%%
    param = rand(342,1);
    param = param/sum(param);
    vals = hist(super_lmnrnd(log(param)', 2000), 1:length(param))';
    p_not_significant = discrete_test(param, vals)

    param_ = param + 0.001*rand(342,1);
    param_ = param_/sum(param_);
    p_significant = discrete_test(param_, vals)

end