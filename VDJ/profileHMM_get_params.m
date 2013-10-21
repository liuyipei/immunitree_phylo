function params = profileHMM_get_params()

    % Chances for having at least i addition, given that we already observed
    % i-1 additions
    ihmmune = {
    '[0]  0.94  [1]  0.85  [2]  0.85  [3]  0.83  [4]  0.8  [5]  0.8  [6]  0.8  [7]  0.8  [8]  0.8  [9]  0.8  [10]  0.8  [11]  0.75  [12]  0.75  [13]  0.75  [14]  0.75  [15]  0.75  [16]  0.75  [17]  0.75  [18]  0.75  [19]  0.75  [20]  0.75  [21]  0.75  [22]  0.75  [23]  0.75  [24]  0.75  [25]  0.75  [26]  0.75  [27]  0.75  [28]  0.75  [29]  0.75  [30]  0.75  [31]  0.75  [32]  0.75  [33]  0.75  [34]  0.75  [35]  0.75  [36]  0.75  [37]  0.75  [38]  0.75  [39]  0.75', ...
    '[0]  0.97  [1]  0.94  [2]  0.92  [3]  0.89  [4]  0.87  [5]  0.85  [6]  0.83  [7]  0.8  [8]  0.8  [9]  0.8  [10]  0.8  [11]  0.8  [12]  0.8  [13]  0.8  [14]  0.8  [15]  0.8  [16]  0.8  [17]  0.8  [18]  0.8  [19]  0.8  [20]  0.8  [21]  0.8  [22]  0.8  [23]  0.8  [24]  0.8  [25]  0.8  [26]  0.8  [27]  0.8  [28]  0.8  [29]  0.8  [30]  0.8  [31]  0.8  [32]  0.8  [33]  0.8  [34]  0.8  [35]  0.8  [36]  0.8  [37]  0.8  [38]  0.8  [39]  0.8'...
    };
    for i=1:2
        ihmmune(i) = textscan(ihmmune{i}, '[%*d] %f');    
    end

    % Chances for having at least i deletions, given that we already observed
    % i-1 deletions.  Values are based on real data.
    % Hence the number of deletion can go from 0 to 21 (22 values).
    ihmmune(3) = {0.57*ones(21,1)}; % V side of VD
    ihmmune(4) = {0.81*ones(21,1)}; % D side of VD
    ihmmune(5) = {0.81*ones(21,1)}; % D side of DJ
    ihmmune(6) = {0.83*ones(21,1)}; % J side of DJ

    for i=1:2
        N_add(i).cdf = [1;exp(cumsum(log(ihmmune{i})))]; % prob (N >= j-1)
        N_add(i).inc = [ihmmune{i}; 0];  % prob (N > j-1 | N >= j-1)
        N_add(i).pdf = N_add(i).cdf .* (1-N_add(i).inc);
        N_add(i).inc2 = [exp(-diff(log(1-N_add(i).cdf))); 1]; % prob (N < j-1 | N <= j-1)
    end

    % N is the number of deletions.  j indexes the vectors
    for i=1:4
        P_del(i).cdf = [1;exp(cumsum(log(ihmmune{i+2})))]; % prob (N >= j-1)
        P_del(i).inc = [ihmmune{i+2}; 0];  % prob (N > i-1 | N >= j-1)
        P_del(i).pdf = P_del(i).cdf .* (1-P_del(i).inc); % prob (N == j-1)
        P_del(i).inc2 = [exp(-diff(log(1-P_del(i).cdf))); 1]; % prob (N < j-1 | N <= j-1)
    end

    % emisison probabilities for the added nucleotide
    N_emit = [0.15 0.35 0.35 0.15];

    params.noise = 0.01;
    params.P_del = P_del;
    params.N_add = N_add;
    params.N_emit = N_emit;
end