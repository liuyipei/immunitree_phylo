% given alignments of reads to V genes, this function filter the V-genes so
% that if every time a read is aligned to x it is also aligned to y
function [redundant no_show] = find_redundant_Vs(Vs, rep)
    %%  Find redundant V ids
    N = size(Vs,1);
    M = size(rep.V,1);
    
    % ix(j,i) = true  iff read j was aligned to rep.V(i)
    ix = false(N,M);
    for i=1:M
        ix(:,i) = sum(Vs == i, 2) > 0;
    end
   
    sum_ix = sum(ix);
    no_show = find(sum_ix == 0);
    players = find(sum_ix > 0);
    
    M_ = length(players);
    
    % I(i,j) = true iff   read aligned to i --> read aligned to j
    I = false(M,M);
    for i = 1:M
        if mod(i,10)==0, i, end
        if sum_ix(i) == 0, continue; end
        for j=i+1:M
            if sum_ix(j) == 0, continue; end
            I(i,j) = all(ix(:,i) <= ix(:,j)); % every time i is listed so is j
            I(j,i) = all(ix(:,j) <= ix(:,i)); % every time j is listed so is i
            if I(i,j) == true && I(j,i) == true, 
                fprintf('%d %d\n', i, j);
                I(i,j) = false;
            end
        end
    end
    % hence, if I(i,j) = true, we can throw away i
    % ====> every non-zero row is thrown away!

    redundant = find(sum(I,2)' > 0 & sum(ix) > 0);
    
    
end


% redundant =
% 
%   Columns 1 through 27
% 
%      6     7    24    26    34    58    65    70    71    74   109   116   135   138   140   146   150   152   153   170   174   176   178   190   204   219   221
% 
%   Columns 28 through 35
% 
%    224   227   253   305   311   315   322   324

% no_show =
% 
%   Columns 1 through 27
% 
%      1     2     3    21    35    38    51    52    53    63    67    68   242   257   258   260   261   262   263   264   265   266   267   269   270   271   272
% 
%   Columns 28 through 54
% 
%    273   274   276   277   278   280   281   282   283   284   285   286   287   288   289   290   291   292   293   296   297   298   300   306   307   309   312
% 
%   Columns 55 through 60
% 
%    313   314   316   317   329   330