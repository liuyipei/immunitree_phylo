function [P H] = multi_hist_with_significance(X,labels, titles, per_patient, interval, normalize, measure1, measure2)
if nargout == 0, figure; end
if ~exist('normalize', 'var'), normalize = 1; end
if ~exist('interval', 'var'), interval = 1; end
if ~exist('per_patient', 'var'), per_patient = false; end
if ~exist('titles', 'var'), titles = cell(size(X,2),1); end
    
T = size(X,2);
P = cell(T,1);
H = cell(T,1);
if isscalar(interval), interval = interval*ones(1,T); end

if per_patient
    M = labels.num+1;
    for i=1:T        
        x = X(:,i);
        tl = titles{i};
        if nargout == 0
            %subplot(T,M,((i-1)*M + (1:(M-1))));
%            hist_with_labels(x, labels, tl, 'clones', interval(i), normalize, 'separate', [T M (i-1)*(M)+1]);
            hist_with_labels(x, labels, tl, 'variants', interval(i), normalize, 'cdfs', [T M (i-1)*(M)+1]);
            subplot(T,M,((i-1)*M + M));
            kstest_with_labels(x, labels, tl, 'per patient');        
            title('2-level random effect test');
        else 
            [P{i} H{i}] = kstest_with_labels(x, labels, tl, 'per patient', measure1, measure2);        
        end
    end
else
    M = 5; %labels.num+1;
    for i=1:T
        x = X(:,i);
        tl = titles{i};
        if nargout == 0
            subplot(T,M,((i-1)*M + (1:(M-1))));
            hist_with_labels(x, labels, tl, 'clones', interval(i), normalize);
            subplot(T,M,((i-1)*M + M));
            kstest_with_labels(x, labels, tl, 'aggregate');        
            title('KS test (aggregated)');
        else 
            [P{i} H{i}] = kstest_with_labels(x, labels, tl, 'aggregate');        
        end
    end
end

end