% tests how likely it is to get the values in x from a uniform[0,1]
% distribution.
% http://en.wikipedia.org/wiki/Fisher's_method
function p = fisher_test(x)
    X2 = -2*sum(log(x(:)));  % higher -> more extreme
    p = 1-chi2cdf(X2,2*length(x)); % P( chi2 >= X2)
end


function test()
%%
    
    junk = zeros(1e4,1);
    junk2 = zeros(1e4,1);    
    for i=1:length(junk)
        x = rand(10,1);        
        junk2(i) = -2*sum(log(x(:)));
        junk(i) = fisher_test(x);
    end
    subplot(2,2,1);
    hist(junk,100)
    title('10,000 fisher test results, null hypothesis');

    subplot(2,2,3);
    bins = 0:0.1:100;
    hist(junk2, bins);
    hold on;
    Z = sum(chi2pdf(bins, 2*10));
    plot(bins, 1e4/Z*chi2pdf(bins, 2*10), 'r')
    title(sprintf('histogram of the fisher statistics over the 10,000 p-values (fisher test score %.2f)', fisher_test(junk)));


    for i=1:length(junk)
        x = 0.95*rand(10,1);
        junk2(i) = -2*sum(log(x(:)));
        junk(i) = fisher_test(x);
    end
    subplot(2,2,2);
    hist(junk,100)
    title('10,000 fisher test results, 0.95*rand(10,1)');

    subplot(2,2,4);    
    bins = 0:0.1:100;
    hist(junk2, bins);
    hold on;
    Z = sum(chi2pdf(bins, 2*10));
    plot(bins, 1e4/Z*chi2pdf(bins, 2*10), 'r');
    title(sprintf('histogram of the fisher statistics over the 10,000 p-values', fisher_test(junk)));

    
    
%%
    
end
