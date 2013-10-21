% returns a histogram based on the data. JV x Source  x patient
% sub_sample = 0: default, no subsample
% sub_sample = 1: for every patient, randomly remove reads from each sample
%                 (except the smallest one) until they are all equal size
% sub_sample = a vector of indexes:  use just those data entries.  
function [big_hist data_] = get_big_hist(data, normalize_per_sample, sub_sample)

if ~exist('sub_sample', 'var'), sub_sample = 0; end

data_ = data;
if ~isscalar(sub_sample) % it's a vector
    % remove naive B cells
    flds = fields(data);
    for f = 1:length(flds)
        data_.(flds{f}) = data_.(flds{f})(sub_sample,:);
    end
end
if sub_sample == 1
    % create a random sub-sampling of the data, so each source has same 
    % sample size (for each individual)
    n = zeros(1,4);
    flds = fields(data);
    for i=1:10 % for each patient   
        for s=1:4
            n(s) = sum(data_.subject_num == (i-1)*4+s);
        end
        nSamples = min(n(n>0));
        for s=1:4
            ix = get_reads(data_,(i-1)*4+s);
            % sub_sample
            iy = randperm(length(ix));
            erase = ix(iy( 1:length(iy)-nSamples));
            for f = 1:length(flds)
                data_.(flds{f})(erase,:) = [];
            end
        end
    end
end

%  Full repertoire per patient and per source
nV = 57; nJ = 6;
big_hist = zeros(nV, nJ, 40);
for k=1:40
    [~, ~, X] = collect_VJs(data_, get_reads(data_,k), 'ihmmune_collapsed', 0);
    if nV > size(X,1), X(nV, nJ) = 0; end    
    if normalize_per_sample
        big_hist(:,:,k) = X./ sum(X(:));  
    else
        big_hist(:,:,k) = X;
    end
end

big_hist = permute(big_hist, [2 1 3]); % JxVx(S*P)
big_hist = reshape(big_hist, [], 4, 10); % (J*V)xSxP


fix_spleen1 = true;
if fix_spleen1
    % fix patient1/spleen to just be the average of the rest
    big_hist(:,2, 1) = mean(big_hist(:,[1 3 4],1), 2);
    if ~normalize_per_sample, big_hist(:,2,1) = round(big_hist(:,2, 1)); end
end




end