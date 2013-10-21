function [tree, sequences, F] = ...
    generate_sequences_from_prior(S, F, rates, Q, pi)

% generate the tree
tree = [0 0 0];
cur_time = 0;
M = 512;  
exprnd_stream = exprnd(1/(rates.birth+rates.death), M,1);

u = 1;
while cur_time < F && sum(tree(:,3)>0) < S
    if u>M, 
        exprnd_stream = [exprnd_stream; exprnd(1/(rates.birth+rates.death), M,1)];
        M = 2*M;
    end
    %next_event = exprnd(1/(rates.birth+rates.death));
    tree(u,3) = min(F, cur_time+exprnd_stream(u));
%    tree(u,3) = cur_time+exprnd_stream(u);
    
    if cur_time > tree(u,2) % if you were not just born
        if rand < rates.death/(rates.death+rates.birth)
            tree(u,3) = -cur_time;
        else
            tree = [tree; u cur_time cur_time];
        end
    end
        
    ix = find(tree(:,3) >= 0);
    if isempty(ix), cur_time = F; break; end
    [cur_time, u] = min(tree(ix,3));
    u = ix(u);
end

tree(:,3) = abs(tree(:,3));
if cur_time<F 
    %fprintf('Finished before time.\n'); end
    F = cur_time;
    tree(:,3) = min(tree(:,3), cur_time);
end

% TODO: check that works for pi coming from synthetic data
sequences = int16(tree_sample_for_phylo_unlog(tree(:,1), size(Q,2)/size(Q,1), Q, pi));
