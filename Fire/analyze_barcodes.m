close all
[text pre post] = textread('/afs/cs/u/joni/scratch/data/lymph/barcodes.txt', '%s %s %s');
text = reshape(text, 73, 6); text = reshape(text(1:40,:)', [], 1); 
pre = reshape(pre, 73, 6); pre = reshape(pre(1:40,:)', [], 1); 
post = reshape(post, 73, 6); post = reshape(post(1:40,:)', [], 1); 



B = pdist(cell2mat(post), 'hamming');
B = length(post{1})*squareform(B);
imagesc(B)
plot_class_lines(1:4:41)
plot_class_lines(1:4:41, true)

C = zeros(10);
for i=1:10,
    for j=1:i
        x = B( (i-1)*4+(1:4), (j-1)*4+(1:4));
        C(i,j) = mean(x(:));
    end
end

figure; imagesc(C);


D = [0     0     0     0     0     0     0     0     0     0;
     4     0     0     0     0     0     0     0     0     0;
     6    20     0     0     0     0     0     0     0     0;
     4     4   102     0     0     0     0     0     0     0;
     8    13     3     6     0     0     0     0     0     0;
     7    12    10     8     6     0     0     0     0     0;
    12    15     9     8     7    12     0     0     0     0;
    13     6     8    10     8     7    22     0     0     0;
    11    18    10     9    13    10    18    13     0     0;
     6     7     5     4     7    11    17    14    10     0];
 
 figure;
 %D = D+D';
 %imagesc(-log(1./(1+D)));
 
 
 %%
 
B2 = pdist(cell2mat(pre), 'hamming');
B2 = length(pre{1})*squareform(B2);
imagesc(B+B2)
plot_class_lines(1:24:241, false, 'm4');
plot_class_lines(1:24:241, true, 'm4');
plot_class_lines(1:6:241, false, 'k2');
plot_class_lines(1:6:241, true, 'k2');

 
%%

