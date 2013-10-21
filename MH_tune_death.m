function [tree codes] = MH_tune_death(tree, F, rates, t)

global greedy

%ix = find(tree(:,1) ~= -1 & tree(:,3)<F);  % all the dead cells (who are not deleted)

can_die = true(size(tree,1),1);  
can_die(t) = false;  % cells that can die are those with no reads
ix = find(tree(:,1) ~= -1 & can_die);

% T <-- the last event in the life of all the sequences in ix
iy = ismember(tree(:,1), ix);  % descendents of nodes in ix
T = sortrows([ix tree(ix,2) ; tree(iy, 1:2)]);
[ix,I,~] = unique(T(:,1),'last');
T = T(I,:);

% remark: notice that we use the rates.birth when we draw the new death time
% That's on purpose!  using this transition probability the acceptance 
% probability is 1, for all cases, *except* when a live cell becomes dead  
% or vise versa.
res = exprnd(1/(rates.birth+rates.death), size(T,1), 1);
res = min(T(:,2)+res, F);
 
% find potentially rejected proposals (iy)
was_alive = tree(ix,3) == F;
now_dead = res<F;
iy = find(was_alive & now_dead); 
iy = [iy ; -find(~was_alive & ~now_dead)];

% find how many cells are alive now (including those with reads)
nAlive = sum(tree(:,3) == F & tree(:,1) ~= -1);

% for any proposal that changes life to death or vise versa:
% if life to death, the acceptance is:
% rates.death/(rates.death+rates.birth)  * (nAlive/(nAlive-1))^N
% (otherwise it's 1/that)

% ix - all cells that are proposed an update
% iy - cells that changed status (alive/dead)
% code(1) = number of cells with same status and different death time
% code(2) = number of cells who were alive and now dead
% code(3) = number of cells who were dead and now alive
% code(4) = number of cells with same status and same death time
% (rejections)
codes = zeros(1,4);
codes(1) = length(ix)-length(iy);  
for k=randperm(length(iy))
    if iy(k) > 0  % a cell died (rare because most live cells have reads)
        acc_prob = rates.death/(rates.death+rates.birth) * (nAlive/(nAlive-1))^length(t);
        if rand<acc_prob  || (greedy == -1 && acc_prob > 0.5)
            codes(2) = codes(2)+1;
            nAlive = nAlive-1;
        else
            codes(4) = codes(4) + 1;
            res(iy(k)) = tree(ix(iy(k)), 3);
        end
    else % a cell born
        acc_prob = (rates.death+rates.birth)/rates.death * (nAlive/(nAlive+1))^length(t);
        if rand<acc_prob || (greedy == -1 && acc_prob > 0.5)
            codes(3) = codes(3)+1;
            nAlive = nAlive+1;
        else
            codes(4) = codes(4) + 1;
            res(-iy(k)) = tree(ix(-iy(k)), 3);
        end
    end      
end

tree(ix,3) = res;


% remark: notice that we use the rates.birth when we draw the new death time
% That's on purpose!  using this transition probability the acceptance 
% probability is 1, because the transition is proportional to the likelihood
% res = exprnd(1/(rates.birth+rates.death), size(T,1), 1);
% res = mod(res, F-T(:,2));
% tree(ix, 3) = T(:,2)+res;
% total_ll_diff = -(rates.birth+rates.death)*(sum(res)-sum(old_res));


end
