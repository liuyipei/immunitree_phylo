function [a b y] = factorize(X)
    % X is a matrix of size(length(a) x length(b)).  
    % goal:  find best a, b such that X - a * b' is minimized.
    
    log_space = false;
    
    [N M] = size(X);    
    
    if log_space
        % option 1: least squares in log-space.
        % that is x_ij = a_i + b_j + epsilon
        A = repmat(eye(N), M,1);
        B = eye(M);
        ix = repmat(1:M, N,1);
        A = [A B(ix(:), :)];

    %     ab = A \ X(:);  
    %     a = ab(1:N);
    %     b = ab(N+1:end);
    %     y = reshape(A*ab, size(X,1), size(X,2));

        ab = A \ log(X(:));  
        a = exp(ab(1:N));
        b = exp(ab(N+1:end));
        y = reshape(exp(A*ab), size(X,1), size(X,2));
    
    else
        % option 2: iteration facorization

        % init a,b
        a = rand(N,1);
        b = rand(M,1);

        nIterations = 20;

        for it = 1:nIterations
            b = (a \ X)';  % a*b' = X
            a = (b \ X')'; % b*a' = X'       
        end
        y = a*b';
    end
end
