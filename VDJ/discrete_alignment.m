function [score editR h]=discrete_alignment(A,R, K)
%function [score editR h]=discrete_alignment(A,R, K)
% Find the probability that sequence A generated read R
% assumption - can't be more than one consecutive insertion/deletion
% answers how A generated R.
% returns the edit string that when applied on R, it aligns to A.
    if nargin == 0, unittest(); return; end
    if nargin<3, K = 5; end
    
    lenA = length(A);
    lenR = length(R);    

    % costs
    eps = 0.001;
    indel = 0.003;
    sigma = 0.02;

    l_del = log(indel);  % letter deleted from A - add insertion to editR
    l_ins = log(indel);  % letter inserted to A - add deletion to editR 
    l_hit = log(1-2*indel); % proceed to next match state
    l_del0 = l_del-l_hit+eps;  % prefix of A deleted - add insertions to editR
    l_ins0 = l_ins-l_hit+eps;  % append to A from left-  add deletions in editR
    
    p_match_to_unobserved_seq = log(0.25);
    p_match_to_unobserved_read = 0;
    p = (1-4*sigma/3)*eye(5);
    p = p + sigma/3;
    p = log(p);
    p(end,:) = p_match_to_unobserved_seq;
    p(:, end) = p_match_to_unobserved_read;
%    q = [log(0.25)*ones(1,4) 0];
    q = log(0.25);

    % main dynamic programming matrix
    fm = -inf*ones(lenA+2,lenR+2);
    fm(2:end,2) = l_del0*(0:lenA)';
    fm(2,2:end) = (q+l_ins0)*(0:lenR);

    % traceback matrix
    bm = zeros(lenA+2,lenR+2);

    %   Going Forward - O(n)
%     for i=1:lenR     % i goes over read (columns of fm)
%         [fm(3:end,i+2) bm(3:end,i+2)] = max(...
%             [fm(1:end-2, i+1) + l_del,...  % letter in A deleted
%              fm(2:end-1, i)   + l_ins + q,...  % insert to A
%              fm(2:end-1, i+1) + l_hit],...  % match
%              [], 2); % max over rows, get a column
%         fm(3:end,i+2) = fm(3:end,i+2) + p(A, R(i));
%     end
    
    for i=3:size(fm,2)     % i goes over read (columns of fm)
        i_minus_k = max(3, i-K);
        i_plus_k = min(size(fm,1), i+K);
        ix  = i_minus_k:i_plus_k;        
        vals = [fm(ix-2, i-1) + l_del,...  % letter in A deleted
                fm(ix-1, i-2) + l_ins + q,...  % insert to A
                fm(ix-1, i-1) + l_hit];  % match
        [fm(ix,i) bm(ix,i)] = max(vals,[], 2); % max over rows, get a column
        fm(ix,i) = fm(ix,i) + p(A(ix-2), R(i-2));
    end

    score = fm(end,end);

    if nargout > 1
        loc = size(fm)';
        steps = [-2  -1  -1;  % sequence
                 -1  -2  -1]; % read
        edit_label = {[1 0], [-1 0], [0]};


        %   Going Backward
%        editR=zeros(1,lenR+lenA);
        editR = [];
        h = [0 0 0];
        id = bm(loc(1), loc(2));

        while (id ~= 0)
            editR = [edit_label{id} editR];
            h(id) = h(id) + 1;
            loc = loc + steps(:,id);
            id = bm(loc(1), loc(2));
        end
        [prefix_count prefix_type] = max(loc-[2 2]');
        editR = [edit_label{prefix_type}(1)*ones(1,prefix_count) editR];
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

[score editR h] = discrete_alignment(A,A);
assert(all(editR == 0));
[score editR h] = discrete_alignment(A,R1);
assert(all(editR == 0));
[score editR h] = discrete_alignment(A,R2);
assert(strcmp('5M-1I-13M-1D-11M', pretty_edit(editR)));
[score editR h] = discrete_alignment(A,R3);
assert(strcmp('5M-1I-15M-1D-9M', pretty_edit(editR)));
[score editR h] = discrete_alignment(A,R4);
assert(strcmp('4D-30M', pretty_edit(editR)));;
[score editR h] = discrete_alignment(A,R5);
assert(strcmp('3I-27M', pretty_edit(editR)));;
fprintf('Test passed!\n');
end


