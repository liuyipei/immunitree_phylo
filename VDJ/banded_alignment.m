function [score editR h]=banded_alignment(A,R,K)
%function [score editR h]=banded_alignment(A,R,K)
    if nargin == 0, unittest(); return; end
    if nargin <3, K = 5; end

    lenA = length(A);
    lenR = length(R);

    if (abs(lenA-lenR) > K)
%        fprintf('warning - length diff > K\n');
        score = -inf;
        editR = []; h = [];
        return;
%        K = abs(lenA-lenR);
    end

    % costs
    eps = 0.001;
    indel = log(0.003);
    sigma = 0.02;

    sx = indel+eps;
    dx = indel;  % insertion in read
    dy = indel;  % deletion in read
    p_match_to_unobserved_seq = log(0.25);
    p_match_to_unobserved_read = 0;
    p = (1-4*sigma/3)*eye(5);
    p = p + sigma/3;
    p = log(p);
    p(end,:) = p_match_to_unobserved_seq;
    p(:, end) = p_match_to_unobserved_read;
    q = [log(0.25)*ones(1,4) 0];
    
    % main dynamic programming matrix
    fm = -inf*ones(lenA+1,lenR+1);
    fm(:,1) = sx*(0:lenA)';
    fm(1,:) = (sx+q(1))*(0:lenR);

    % traceback matrix
    bm = zeros(lenA+1,lenR+1);

    %   Going Forward - O(n)
    for i=2:size(fm,1)     % i goes over sequence
        i_minus_k = max(2, i-K);
        i_plus_k = min(size(fm,2), i+K);
        for j=i_minus_k:i_plus_k % j goes over read

            % coming from the top
            vals = [fm(i-1,j) + dx;             % from top, insertion in read
                    fm(i,j-1) + dy + q(R(j-1));   % from left, deletion in read
                    fm(i-1,j-1)+p(A(i-1),R(j-1));  % diag
                    ];
            [fm(i,j) bm(i,j)] = max(vals);
        end
    end


    score = fm(end,end);

    if nargout > 1
        loc = size(fm)';
        steps = [-1  0 -1;
                 0  -1 -1];
        edit_label = [1 -1 0];


        %   Going Backward
        editR=zeros(1,lenR+lenA);
        e_ix = lenR + lenA;
        h = [0 0 0];

        id = bm(loc(1), loc(2));

        while (id ~= 0)    
            editR(e_ix) = edit_label(id); e_ix = e_ix -1;
            h(id) = h(id) + 1;
            loc = loc + steps(:,id);
            id = bm(loc(1), loc(2));
        end
        [prefix_count prefix_type] = max(loc-[1 1]');
        editR = [edit_label(prefix_type)*ones(1,prefix_count) editR(e_ix+1:end);];
    end
end


function unittest()
%%
A =  [1 2 3 2 2 2 4 4 3 1 1 2 3 2 1 4 1 2 2   3 2 1 4 2 2 4 2 1 1 1];
R1 = [1 2 3 2 2 2 4 4 3 1 2 2 3 2 1 4 1 2 2   2 2 1 4 2 2 4 2 1 1 1];
R2 = [1 2 3 2 2   4 4 3 1 1 2 3 2 1 4 1 2 2 2 3 2 1 4 2 2 4 2 1 1 1];
R3 = [1 2 3 2 2   4 4 3 1 2 2 3 2 1 4 1 2 2 2 2 2 1 4 2 2 4 2 1 1 1];
R4 = [3 1 1 1 1 2 3 2 2 2 4 4 3 1 1 2 3 2 1 4 1 2 2   3 2 1 4 2 2 4 2 1 1 1];
R5 = [2 2 2 4 4 3 1 1 2 3 2 1 4 1 2 2   3 2 1 4 2 2 4 2 1 1 1];

K = 4;
[score editR] = banded_alignment(A,A, K);
assert(all(editR == 0));
[score editR] = banded_alignment(A,R1, K);
assert(all(editR == 0));
[score editR] = banded_alignment(A,R2, K);
assert(strcmp('5M-1I-13M-1D-11M', pretty_edit(editR)) || strcmp('3M-1I-13M-1D-13M', pretty_edit(editR)));
[score editR] = banded_alignment(A,R3, K);
assert(strcmp('5M-1I-15M-1D-9M', pretty_edit(editR)));
[score editR] = banded_alignment(A,R4, K);
assert(strcmp('4D-30M', pretty_edit(editR)));;
[score editR] = banded_alignment(A,R5, K);
assert(strcmp('3I-27M', pretty_edit(editR)));;
fprintf('Test passed!\n');

%%
[score editR] = banded_alignment(A,A, K)
[score editR] = banded_alignment(A,R1, K)
[score editR] = banded_alignment(A,R2, K)
[score editR] = banded_alignment(A,R3, K)
[score editR] = banded_alignment(A,R4, K)
[score editR] = banded_alignment(A,R5, K)

end


