function [nMuts mut mut_str] = annotate_mutations_on_tree(...
    T, sequences, gap_is_mutation, germline_count_end_indx, do_pairwise_mutstr)

if ~exist('gap_is_mutation', 'var'), gap_is_mutation = false; end
if ~exist('germline_count_end_indx', 'var'), germline_count_end_indx = Inf; end
if max(sequences(:)) > 5, sequences = codons2seqs(sequences); end
if ~exist('do_pairwise_mutstr', 'var'), do_pairwise_mutstr = false; end

M = size(sequences,1);
if isequal(size(T), [1 size(sequences,2)])
    sequences = [T; sequences];
    T = [0; ones(M,1)];
    M = M+1;
end
nMuts = zeros(M, 1);
mut_str = cell(M, 1);
mut = -2*ones(size(sequences));
map = 'ACGT-';
for t=1:M
    mut_str{t} = '';
    if T(t,1) == 0, continue; end % germline not annotated
    
    p = T(t, 1);
    x1 = sequences(t,:);
    b1 = sequences(p,:); % parent_seq
    
    
    % pairwise alignment for just these two -- so we can avoid cases
    % like 100A->_  100_->A in the same string
    
    % use the global scheme for nMuts, mut
    % use the pairwise scheme for mutation strings
    ignore_lead_gaps = true;
    if ~do_pairwise_mutstr
        x = x1; b = b1;
    else
        glocalvalue = gap_is_mutation & ignore_lead_gaps;
        [x2, b2, aligned2] = prettyalign2seqs(x1,b1,glocalvalue);
        x = x2; b = b2;
    end
    
    if gap_is_mutation        
        
        if ignore_lead_gaps
            % 02/11/2012 -- Yi -- primarily thinking about
            % ignoring description for lead gaps for nodes connected to germline
            mut_loci = find(b ~= x);
            first_non_gap_pos = find(b ~= 5 & x ~= 5, 1, 'first');
            if length(first_non_gap_pos) < 1
                first_non_gap_pos = 0;
            end
            mut_loci = mut_loci(mut_loci >= first_non_gap_pos);
        else
            mut_loci = find(b ~= x);
        end
    else
        mut_loci = find(b ~= x & b ~= 5);
    end
    
    if (p == 1)
        mut_loci = mut_loci(mut_loci <= germline_count_end_indx);  % defaults to infinity, except when the parent is 1
    end
    mut_val = x(mut_loci);
    parent_val = b(mut_loci);
    
    if ~do_pairwise_mutstr
        mut(t, mut_loci) = mut_val;
    else
        mut = [];
    end
    nMuts(t) = length(mut_loci);
    %if nMuts(t) == 0, nMuts(t) = -1; end % not to ruin the graph
    parent_seq_index = cumsum((b ~= 5));
    mut_str{t} = sprintf('%d%c>%c|', ...
        [parent_seq_index(mut_loci); map(parent_val);map(mut_val)]);
    if ~isempty(mut_str{t})
        mut_str{t} = mut_str{t}(1:end-1); % remove redundant '|'  delimiter
    end

end % end for each node

end

function test()
%%
[nMut mut mut_str] = annotate_mutations_on_tree(tree_, codons2seqs(sequences_));

end

function [seq1, seq2, aligned] = prettyalign2seqs(seq1, seq2, glocalvalue)
    % convert numbers into letters, do alignment, and convert back into numbers
    map_from_int = 'ACGT-';
    map_to_int = [];
    map_to_int('ACGTNRacgtn-') = [1:5 5 1:5 5];
    seq1_ = map_from_int(non5(seq1));
    seq2_ = map_from_int(non5(seq2));
    [~, aligned] = nwalign(seq1_, seq2_, 'glocal', glocalvalue);
    seq1 = map_to_int(aligned(1,:));
    seq2 = map_to_int(aligned(3,:)); % second row is colons, pipes, spaces
end
function non5 = non5(x)
    non5 = x(x~=5);
end