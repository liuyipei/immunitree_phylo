function big_hist_permuted = permute_big_hist(big_hist, permutation_code)

big_hist = reshape(big_hist, 57*6, 4, 10);
big_hist_permuted = big_hist;

switch permutation_code
    case 2 % switch sources
        for k=1:10
            ix = randperm(4);
            big_hist_permuted(:,:,k) = big_hist_permuted(:,ix,k);
        end
  

    case 3 % switch non-blood sources
        for k=1:10  
            ix = 1+randperm(3);
            big_hist_permuted(:,2:4,k) = big_hist_permuted(:,ix,k);
        end
        
    case 4  % switch patients
        for k=1:4  
            ix = randperm(10);
            big_hist_permuted(:,k,:) = big_hist_permuted(:,k,ix);
        end        
    
    case 5 % switch VJ combinations
        ix = randperm(size(big_hist,1));
        big_hist_permuted = big_hist_permuted(ix,:,:);

end
    
    

end