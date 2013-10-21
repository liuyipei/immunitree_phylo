function [rmse S P big_hist_ delta] = factorize_big_hist(big_hist, rep, vj_factorization, mean_factorization)    
    if ~exist('vj_factorization', 'var'), vj_factorization = false; end

    nV = length(rep.V);
    nJ = length(rep.J);
    nS = size(big_hist,2);
    nP = size(big_hist,3);
    

    big_hist_baseline = mean(reshape(big_hist, size(big_hist,1), []), 2);
    big_hist_baseline = big_hist_baseline(:, ones(1,4), ones(1,10));
    
    if mean_factorization == 1
        big_hist = big_hist - big_hist_baseline;
    elseif mean_factorization == 2        
        big_hist = big_hist ./ big_hist_baseline;
        big_hist(isnan(big_hist)) = 0;
    end
    
    big_hist_ = big_hist;

    S = zeros(size(big_hist_,2), size(big_hist_,1));
    P = zeros(size(big_hist_,3), size(big_hist_,1));
    for i=1:size(big_hist_,1)
        x = squeeze(big_hist_(i, :, :));
        x = x+1e-6*randn(size(x));
        [S(:,i) P(:,i) big_hist_(i, :, :)] = factorize(x);      
    end

    % further factorization of S and p
    
    if vj_factorization
        S = reshape(S, [], nJ, nV);
        P = reshape(P, [], nJ, nV);

        SV = zeros(nV, nS);
        SJ = zeros(nJ, nS);
        PV = zeros(nV, nP);
        PJ = zeros(nJ, nP);
    
        % factorize each source repertoire S(i,:,:)
        for i=1:size(S,1)
            x = reshape(S(i,:,:), nJ, nV);
            x = x+1e-6*randn(size(x));
            [SJ(:,i) SV(:,i) S(i,:,:)] = factorize(x);      
        end
        
        % factorize each patient repertoire P(i,:,:)
        for i=1:size(P,1)
            x = reshape(P(i,:,:), nJ, nV);
            x = x+1e-6*randn(size(x));
            [PJ(:,i) PV(:,i) P(i,:,:)] = factorize(x);      
        end
        
        % reconstruct S,P and big_hist_ using small factors
        S = reshape(S, size(big_hist_,2), []);
        P = reshape(P, size(big_hist_,3), []);    
        for i=1:size(big_hist_,1)
            big_hist_(i, :, :) = S(:,i)*P(:,i)';
        end
    end
    
%     rmse = sqrt(mean( (big_hist(:)-big_hist_(:)).^2))
    
    if vj_factorization %even_more_factorization
        S = reshape(S, [], nJ, nV);
        P = reshape(P, [], nJ, nV);
        big_hist  = reshape(big_hist,  nJ, nV, nS, nP);
        big_hist_ = reshape(big_hist_, nJ, nV, nS, nP);
        
        % coordinate optimization
        
        for it = 1:10
        % optimize pV:
        % for each patient p, gene V
        for p=1:nP
            for i=1:nV
                % collect all entries in the big_hist relevant to (p,v) -> x
                % and the corresponding ones in big_hist_ -> y
                % divide y by current pV(p,v)        
                % new pV(p,v) = x'*y / norm(y)^2
                old_val = PV(i,p);
                x = big_hist(:,i,:,p);
                y = big_hist_(:,i,:,p); 
                if ~any(y(:)), continue; end
                val = old_val * ( x(:)'*y(:)) / ( y(:)'*y(:));
                big_hist_(:,i,:,p) = (val/old_val) * y;
                PV(i,p) = val;
            end

            for i=1:nJ
                old_val = PJ(i,p);
                x = big_hist(i,:, :,p);
                y = big_hist_(i,:,:,p); 
                if ~any(y(:)), continue; end
                val = old_val * ( x(:)'*y(:)) / ( y(:)'*y(:));
                big_hist_(i,:,:,p) = (val/old_val) * y;
                PJ(i,p) = val;
            end
            P(p,:,:) = PJ(:,p)*PV(:,p)';
    
        end

        
        for s=1:nS
            for i=1:nV
                old_val = SV(i,s);
                x = big_hist(:,i,s,:);
                y = big_hist_(:,i,s,:); 
                if ~any(y(:)), continue; end
                val = old_val * ( x(:)'*y(:)) / ( y(:)'*y(:));
                if isnan(val) || (old_val == 0), keyboard; end
                big_hist_(:,i,s,:) = (val/old_val) * y;
                SV(i,s) = val;
            end

            for i=1:nJ
                old_val = SJ(i,s);
                x = big_hist(i,:,s,:);
                y = big_hist_(i,:,s,:); 
                if ~any(y(:)), continue; end
                val = old_val * ( x(:)'*y(:)) / ( y(:)'*y(:));
                if isnan(val) || (old_val == 0), keyboard; end
                big_hist_(i,:,s,:) = (val/old_val) * y;
                SJ(i,s) = val;
            end
            S(s,:,:) = SJ(:,s)*SV(:,s)';            
        end
%         fprintf('%d: %.g\n', it,sqrt(mean( (big_hist(:)-big_hist_(:)).^2))-rmse);
%         rmse = sqrt(mean( (big_hist(:)-big_hist_(:)).^2));
        end
        
        big_hist  = reshape(big_hist,  [], nS, nP);
        big_hist_ = reshape(big_hist_, [], nS, nP);
%         sanity = big_hist_;
%         for i=1:size(big_hist_,1)
%             big_hist_(i, :, :) = S(:,i)*P(:,i)';
%         end
%         max(abs(big_hist_(:)-sanity(:)))
  
        S = reshape(S, size(big_hist_,2), []);
        P = reshape(P, size(big_hist_,3), []);    
     
    end
    
    if mean_factorization == 1
        big_hist  = big_hist  + big_hist_baseline;
        big_hist_ = big_hist_ + big_hist_baseline;
    elseif mean_factorization == 2        
        big_hist  = big_hist  .* big_hist_baseline;
        big_hist_ = big_hist_ .* big_hist_baseline;
    end    

    Z  = plot_patient_repertoire(big_hist, rep.V, rep.J);
    Z_ = plot_patient_repertoire(big_hist_, rep.V, rep.J);
    delta = Z-Z_;
    rmse = sqrt(mean((delta(:)).^2));
    

end


function [rmse a b big_hist_ delta] = factorize_big_hist2(big_hist, rep)    




end
