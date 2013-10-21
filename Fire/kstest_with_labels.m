function [P H] = kstest_with_labels(x, labels, tl, group, measure1, measure2)

if ~exist('measure2', 'var'), measure2 = 'random_effects'; end    
if ~exist('measure1', 'var'), measure1 = 'random_effects'; end
if ~exist('group', 'var'), group = 'aggregate'; end
if ~exist('tl', 'var'), tl = ''; end

not_nan = ~isnan(x);

if strcmp(group, 'aggregate')
    P = -ones(labels.num, labels.num);
    for l=1:labels.num
        iy = labels.l == l;
        for l_=1:labels.num
            iy_ = labels.l == l_;
            if sum(iy) == 0 || sum(iy_) == 0, continue; end
            if all(isnan(x(iy))) || all(isnan(x(iy_))), continue; end
            [~, P(l,l_)] = kstest2(x(iy),x(iy_), 0.05, 'larger');
        end
    end  
    H = [];

    if nargout == 0
        %figure; 
%        imagesc(-log(P), [0 8]);
        imagesc(P, [0 0.2]);
        colormap summer
        set(gca, 'YTick', 1:labels.num);
        set(gca, 'XTick', 1:labels.num);
        set(gca, 'XTickLabel', cellfun(@(x) x(1), labels.names)');
        set(gca, 'YTickLabel', labels.names);
        title(sprintf('one-sided KS test over empirical distributions of ''%s'', \nfor every pair of ''%s'', aggregated over patients', ...
            tl, labels.str));
    end
    
elseif strcmp(group, 'per patient')
   %% For each patient calcualte SK stats how the distributions are different
    L = max(labels.patients);
    H = nan(labels.num, labels.num, L);
    W = zeros(labels.num, labels.num, L);
    P = ones(labels.num, labels.num);
    for l=1:labels.num
        iy = labels.l == l;
        st = 1;
        if strcmp(measure1, 'random_effects'), st = l+1; end
        for l_=st:labels.num
            iy_ = labels.l == l_;
            for patient=1:L
                iz = labels.patients == patient; 
                w1 = sum(~isnan(x(iz & iy & not_nan)));
                w2 = sum(~isnan(x(iz & iy_ & not_nan)));
                if w1 == 0 || w2 == 0, continue; end
                
                % base measure
                if strcmp(measure1, 'random_effects') % TODO: remove NaNs.
                    [H(l,l_, patient) W(l,l_,patient)] = study_effect(x(iz & iy& not_nan), x(iz & iy_& not_nan));
                elseif strcmp(measure1, 'ttest')
                    [~, H(l,l_, patient)] = ttest2(x(iz & iy& not_nan), x(iz & iy_), 0.05, 'left');                    
                    W(l,l_,patient) = min(w1,w2);
                elseif strcmp(measure1, 'kstest')
                    [~, H(l,l_, patient)] = kstest2(x(iz & iy& not_nan), x(iz & iy_& not_nan), 0.05, 'larger');                    
                    W(l,l_,patient) = min(w1,w2);
                end                                    
            end 
            iw = ~isnan(H(l,l_,:));
            
            % meta measure
            if strcmp(measure2, 'random_effects')
                assert(strcmp(measure1, 'random_effects'));
                P(l,l_) = random_effect(H(l,l_,iw), W(l,l_,iw));
            elseif strcmp(measure2, 'fisher')
                P(l,l_) = fisher_test(H(l,l_,iw));
            elseif strcmp(measure2, 'mytest')
                P(l,l_) = mytest(H(l,l_,iw));
            end
        end
    end   
%    P = Ztest(H, W);
    if nargout == 0        
        %figure; 
%        imagesc(-log(P), [0 8]);
        imagesc(P, [0 0.2]);
        colormap summer
        set(gca, 'YTick', 1:labels.num);
        set(gca, 'XTick', 1:labels.num);
        set(gca, 'XTickLabel', cellfun(@(x) x(1), labels.names)');
        set(gca, 'YTickLabel', labels.names);
%         title(sprintf('fisher test over p-valued from per-patient one-sided KS test\nover empirical distributions of ''%s'', for every pair of ''%s''', ...
%             tl, labels.str));         
    end
else
    fprintf('Error in ''group'' parameter!\n');    
end


end

function P = mytest(H)
    thresh = 0.05;
    N = length(H);
    K = sum(H<thresh); % successful
    P = 0;
    for k=K:N
        P = P + nchoosek(N,k) * thresh^k * (1-thresh)^(N-k);
    end
    P = K;
end
% function P = mytest(H)
%     thresh = 0.05;
%     N = sum(H>-1, 3);
%     K = sum(H >-1 & H<thresh, 3); % successful
%     P = zeros(size(N));
%     for i=1:length(K)
%         for j=1:length(K)
%             for k=K(i,j):N(i,j)
%                 P(i,j) = P(i,j) + nchoosek(N(i,j),k) * thresh^k * (1-thresh)^(N(i,j)-k);
%             end
%         end
%     end
% end

% estimate T = E[avg(x)-avg(y)] and v = Var[avg(x)-avg(y)]
function [T v] = study_effect(x, y)
    assert(~isempty(x) && ~isempty(y));
    assert(all(~isnan(x)) && all(~isnan(y)) );
    
    N = length(x);
    M = length(y);
    
    combined_std = sqrt(  ( (N-1)*var(x) + (M-1)*var(y) ) / (N+M-2)  );
    
    % Cohen's d
    T = (mean(y) - mean(x))/combined_std;
    v = (1/N + 1/M) + T^2 / (2*(N+M));
    
    % Hedges' g
    df = N+M-2;
    J = 1-3/(4*df-1);
    T = J * T;
    v = J^2 * v;

    
%     % Difference of means
%     T = mean(y) - mean(x);
%     v = var(x)/length(x) + var(y)/length(y);
end
