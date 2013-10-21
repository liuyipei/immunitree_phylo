function [aligned_prefix trimmed st sc] = trim_and_fix_reads(X, V, direction, v_cut)
%function [aligned_prefix trimmed st] = trim_and_fix_reads(X, V, direction, v_cut)
% Input: 
% X  a set of reads X, 
% V  a prefix or a suffix to look for 
% direction   is it a prefix (default) or a suffix.
% v_cut   from which point to trim.
% Output: 
% aligned_prefix  the part of read, always with length v_cut-1, that was
%                 aligned to V(1:v_cut-1) (after alignment)
% trimmed         the rest of the read, unaligned
% st              to which letter of V does the first letter of the trimmed 
%                 read corresponds

%    if ~exist('v_cut', 'var'), v_cut = 1; end
    if nargin < 4, v_cut = 1; end
    
    N = length(X);
    st = zeros(2, N);
    sc = zeros(N, 1);
    en = zeros(2, N);
%    trimmed = cell(N,1);
    aligned_prefix = cell(N,1);
    trimmed = cell(N,1);
    
% preperation    
    eps = 0.001;
    indel = 0.003;
    sigma = 0.01;
    
    p_match_to_unobserved_seq = log(0.25);
    p_match_to_unobserved_read = 0;
    p = (1-4*sigma/3)*eye(5);
    p = p + sigma/3;
    p = log(p);
    p(end,:) = p_match_to_unobserved_seq;
    p(:, end) = p_match_to_unobserved_read;
    q = log(0.25);
    penalty = -(log(indel)+q);
    corr = 3;
    map('ACGTNRacgtn-') = [1:5 5 1:5 5];
    V = map(V);
    if direction == 1
        V = fliplr(V);
    end
    
    for i=1:N

        if mod(i,1000) == 0, fprintf('%d ', i); end
        if X{i}(1) > 5, X{i} = map(X{i}); end
        if direction == 1, X{i} = fliplr(X{i}); end
        
        % align V to X{i} using SW.
        [sc(i) al st(:,i)] = swalign(V, X{i}, 'alphabet', 'nt', 'ScoringMatrix', p+corr, ...
            'GapOpen', penalty, 'ExtendGap', penalty-eps-corr);
    
        % divide reads to aligned part and unaligned part

        % aligne prefix of read to V
        non_gaps_in_V = al(1,:) ~= '-';
        aligned = al(3,non_gaps_in_V);  % erase gaps in V
        
        
        % add gaps to R so it will match missing parts of V
        left_pad_v = st(1,i)-1;
        right_pad_v = v_cut-1-(length(aligned) + left_pad_v);
        % aligned part:  read prefix is aligned to V(1:v_cut-1)
        aligned = ['-'*ones(1,left_pad_v) aligned '-'*ones(1,right_pad_v)];  
        aligned_prefix{i} = map(aligned(1:v_cut-1));

        % the rest of the read, unaligned
        % where in the alignment begins the rest of the read?
        v_cut_with_gaps = find(cumsum(non_gaps_in_V) == v_cut-left_pad_v, 1);
        if isempty(v_cut_with_gaps)
            % either the read does not cover the cut point
            if v_cut-left_pad_v >= 0
                v_cut_with_gaps = size(al, 2)+1; 
            else % or the read starts after the cut point
                v_cut_with_gaps = 1;
            end
        end
        al_prefix = al(:,1:v_cut_with_gaps-1);

        % how many read letters were emitted until the V cut point?
        read_letters_observed_in_prefix = sum(al_prefix(3,:) ~= '-');
        % rest of read:
        trimmed{i} = X{i}( (st(2,i)-1+read_letters_observed_in_prefix+1):end);

%       edit string stuff        
%         editR = zeros(1, size(al_prefix,2));
%         editR(al_prefix(3, :) == '-') = 1;
%         editR(al_prefix(1, :) == '-') = -1;        
        
        % if the last letter in the prefix is 5 and 
        % and nothing more is aligned to V, replace insertion with mismatch
%         if ~isempty(aligned_prefix{i}) && aligned_prefix{i}(end) == 5 && ...
%                 length(aligned) == v_cut-1 && ~isempty(trimmed{i})
%             aligned_prefix{i}(end) = trimmed{i}(1);
%             trimmed{i} = trimmed{i}(2:end);
%          end
            
    end
    
    for i=1:N
        if direction == 1, 
            aligned_prefix{i} = fliplr(aligned_prefix{i});
            trimmed{i} = fliplr(trimmed{i});
        end
    end
end

function test()
%%
V =  'ACACAGTG';
X = {'ACACCAC'};
[A T st sc] =trim_and_fix_reads(X, V, 0, 6);
A{1}, T{1}
st
% in this case we should not end with an insertion
assert(A{1}(end) ~= 5);
fprintf('Passed.\n');
%%
V =  'ACACAGTG';
X = {'ACACGTG'};
[A T st sc] =trim_and_fix_reads(X, V, 0, 6);
A{1}, T{1}
st
% in this case we should end with an insertion
assert(A{1}(end) == 5);
fprintf('Passed.\n');
end
