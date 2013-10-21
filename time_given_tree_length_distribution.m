function [Mu Sigma] = time_given_tree_length_distribution(N, N_)

    dN = N_-N+1;
    H = 1./[N:N_]';
    Mu = sum(H)/dN;
    Sigma = sqrt(2/(dN^2+dN) *cumsum(H)'*H - (sum(H)/dN)^2);

end

function test()
%%
N = 20;
N_ = 100;
dN = N_-N+1;

M = 10000;
psi = drchrnd(ones(1,dN), M);

H = 1./[N:N_]';
T = psi*H; 

[Mu Sigma] = time_given_tree_length_distribution(N, N_)
mean(T) - Mu
std(T) - Sigma

figure;
x = 0:0.0001:max(T);
y = M*diff(normcdf(x, sum(H)/dN, std(T)));
hist(T, x);
hold on;
plot(x(1:end-1), y, 'r');


end


% mu = 10;  sigma = 5;
% N = 1e5;
% y = zeros(1,N);
% for i=1:N
%     x = mu+sigma*randn;
%     y(i) = exprnd(x);
% end
% hist(y, 0:0.01:100);

