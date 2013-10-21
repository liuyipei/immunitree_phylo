function [uCA,ndx,pos] = uniqueRowsCA_filler(iCA,TREAT_NAN_EQUAL,FIRST_LAST)
% uses a filler value to pad each row so that you can run uniqueRowsCA on
% it
    filler = -1000;
    maxLen = max(cellfun(@length,iCA));
    iCA2 = cellfun(@(x)[x filler*ones(1,maxLen-length(x))], iCA, 'UniformOutput', false);
    [uCA2,ndx,pos] = uniqueRowsCA(iCA2,TREAT_NAN_EQUAL,FIRST_LAST);    
    uCA = cellfun(@(x)x(x~=filler), uCA2, 'UniformOutput', false);
end
