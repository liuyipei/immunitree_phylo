% input:  estimations of the mean and varaiances over a random variable S,
% obtain from different studies.  T(i) and v(i) are the mean and variance
% for S obtained from study i.  The output is the p-value for the 
% alternative assumption that S ~= 0.  
% fixed = true:  assume one S
% fixed = false: assume for each study the mean is transformed via a a 
% random study specific effect: S+R_i
function [p z T_star V_star] = random_effect(T, v, fixed)
    T = T(:); v = v(:);
    if nargin == 3 && fixed
        [p z T_star V_star] = fixed_effect(T,v);
        return;
    end
        
    k = length(T); 
    assert(length(v) == k);
    w = 1./v;
    T_dot = w'*T / sum(w);
    Q = w'* ((T - T_dot).^2);
    C = sum(w) - (sum(w.^2)/sum(w));
    gamma2 = max(0, (Q-k+1)/C);
    v_ = v +gamma2;
    [p z T_star V_star] = fixed_effect(T,v_);
end

function [p z T_dot v_dot] = fixed_effect(T, v)
    k = length(T); 
    assert(length(v) == k);
    w = 1./v;
    T_dot = w'*T / sum(w);
    v_dot = 1/sum(w);
    z = T_dot / sqrt(v_dot);
    p = 2*(1-normcdf(abs(z)));
end

function test()
%% testing using the numbers from 
% http://www.meta-analysis.com/downloads/Meta-analysis_fixed_effect_vs_random_effects_sv.pdf

T = [0.1  0.3  0.35 0.65  0.45 0.15]';
v = [0.03 0.03 0.05 0.010 0.05 0.02]';
[p,z,T_star,v_star] = random_effect(T, v, true)
assert(abs(z-6.3563)<1e-4);
[p,z,T_star,v_star] = random_effect(T, v, false)
assert(abs(z-3.2247)<1e-4);
fprintf('Test passed!\n');

end