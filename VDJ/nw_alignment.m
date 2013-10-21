function [score editR h]=nw_alignment(A,R, indel, noise)
%function [score editR h]=banded_alignment(A,R,K)
    if nargin == 0, unittest(); return; end

    lenA = length(A);
    lenR = length(R);

    % costs
    eps = 0.001;      
    
    p_match_to_unobserved_seq = log(0.25);
    p_match_to_unobserved_read = 0;
    p = (1-4*noise/3)*eye(5);
    p = p + noise/3;
    p = log(p);
    p(end,:) = p_match_to_unobserved_seq;
    p(:, end) = p_match_to_unobserved_read;
    q = log(0.25);
    penalty = -(log(indel)+q);

    [score al] = nwalign(A, R, 'alphabet', 'NT', 'ScoringMatrix', p, ...
        'GapOpen', penalty, 'ExtendGap', penalty-eps);
    
    editR = zeros(1, size(al,2), 'int8');
    editR(al(3, :) == '-') = 1;
    editR(al(1, :) == '-') = -1;
end    
    

function unittest()
%%
A =  [1 2 3 2 2 2 4 4 3 1 1 2 3 2 1 4 1 2 2   3 2 1 4 2 2 4 2 1 1 1];
R1 = [1 2 3 2 2 2 4 4 3 1 2 2 3 2 1 4 1 2 2   2 2 1 4 2 2 4 2 1 1 1];
R2 = [1 2 3 2 2   4 4 3 1 1 2 3 2 1 4 1 2 2 2 3 2 1 4 2 2 4 2 1 1 1];
R3 = [1 2 3 2 2   4 4 3 1 2 2 3 2 1 4 1 2 2 2 2 2 1 4 2 2 4 2 1 1 1];
R4 = [3 1 1 1 1 2 3 2 2 2 4 4 3 1 1 2 3 2 1 4 1 2 2   3 2 1 4 2 2 4 2 1 1 1];
R5 = [2 2 2 4 4 3 1 1 2 3 2 1 4 1 2 2   3 2 1 4 2 2 4 2 1 1 1];

indel = 0.01;
noise = 0.01;

[score editR] = nw_alignment(A,A, indel, noise);
assert(all(editR == 0));
[score editR] = nw_alignment(A,R1, indel, noise);
assert(all(editR == 0));
[score editR] = nw_alignment(A,R2, indel, noise);
assert(sum(editR~=0) == 2 && sum(editR) == 0);
[score editR] = nw_alignment(A,R3, indel, noise);
assert(sum(editR~=0) == 2 && sum(editR) == 0);
[score editR] = nw_alignment(A,R4, indel, noise);
assert(strcmp('4D-30M', pretty_edit(editR)));
[score editR] = nw_alignment(A,R5, indel, noise);
assert(strcmp('3I-27M', pretty_edit(editR)));
fprintf('Test passed!\n');

end