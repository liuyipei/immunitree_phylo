A = [2 54 6 71;
    4 96 6 77;
    4 96 6 83; 
    4 96 6 89; 
    8 77 6 104; 
    16 77 6 74;
    16 77 6 83;
    28 187 6 80; 
    32 54 12 77; 
    30 77 6 95; 
    32 77 6 89; 
    26 217 6 89];

fprintf('cp ');
for i=1:size(A,1)
    fprintf('sample_%d_V_%d_J_%d_len_%d_clone_1.jpg ', A(i,1), A(i,2), A(i,3), A(i,4));
end
fprintf('~/www/research/VDJ/Boyd_Fire/\n');

for i=1:size(A,1)
    fprintf(' %2d\t%2d\t%d\t%3d\n', A(i,1), A(i,2), A(i,3), A(i,4));
    fprintf('http://ai.stanford.edu/~joni/research/VDJ/Boyd_Fire/sample_%d_V_%d_J_%d_len_%d_clone_1.jpg\n', A(i,1), A(i,2), A(i,3), A(i,4));    
end
   
%%

N= length(data.reads);
vdj_hist = hist(data.ihmmune(:,[1 4 5 6 7]), 0:max(data.ihmmune(:)));
vdj_hist(1,2:end) = 0; 
subplot(4,1,1:2);
bar(vdj_hist, 'stacked');
title('V alignment popularity');
fprintf('%.1f%% of all reads do not have V alignment\n', 100*sum(data.ihmmune(:,1) == 0)/N);

%% 
M = max(data.ihmmune);
I = false(M,M);
for i=1:M
    ix = find(sum(data.ihmmune(:,[1 4:7]) == i,2));
    if isempty(ix), continue; end
    i
    for j=1:M
        junk = sum(data.ihmmune(ix,[1 4:7]) == j,2);
        if all(junk==1), I(i,j) = true; end
    end
end    
I(sub2ind([M M], 1:M, 1:M)) = false;

%%
erase = [];
[J K] = find(I_);
for i=1:length(J)
    if I(K(i), J(i)) == false,
        erase = [erase J(i)];
    end
end

%    3     6     7    24    26    34    58    63    64    70    71    74   109
% 
%    116   135   138   140   144   150   153   170   174   176   178   190   204
% 
%    219   221   226   269   278   309   315   322