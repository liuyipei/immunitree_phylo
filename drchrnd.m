function r = drchrnd(a,n)
    p = length(a);
    r = randg(a(ones(1,n),:));
    sum_r = sum(r,2);
    r = r ./ sum_r(:,ones(1,p)); % normalize
end