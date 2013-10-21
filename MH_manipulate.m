function [tree total_ll_diff codes] = MH_manipulate(tree, sequences, Q, rates, t)

global greedy

total_ll_diff = 0;
codes = zeros(1,5);

B = size(Q,1); 

%%%
M=size(tree,1);
% Using DG is slower
%DG = sparse(1+tree(:,1), 2:M+1, true, M+1, M+1);
%DG = DG(2:end, 2:end);
children = cell(M,1);
for i=1:M
%    children{i} = find(DG(i,:))';
    children{i} = find(tree(:,1) == i);
end
%%%

evidence = nodes_with_evidence(tree,t);
accepted = true;

%for v=randperm(size(tree,1))
for i=1:sum(evidence)
    if accepted, pool = find(evidence); end
    v = pool(ceil(rand*length(pool)));   
    if tree(v,1) <= 0, continue; end  % doesn't handle root or deleted nodes.
    code = [v 0 0];
    [tree_ transition_prob code] = manipulate(tree, code, children);
    if code(1) == 0, codes(5) = codes(5)+1; continue; end
    assert(code(1) ~= 0);    
    assert(code(2) ~= 0);

    bias = 1;
    ll_diff = -(rates.death+rates.birth)*(tree(v,2)-tree_(v,2));
    a = tree(v,1);
    b = tree_(v,1);
    if b ~= a
        ix = find(sequences(b,:) ~= sequences(a,:));
        if ~isempty(ix)
            ll_diff = ll_diff + ...
                -sum(log(Q(mysub2ind(B, sequences(v,ix), sequences(a,ix), ix)))) ...
                +sum(log(Q(mysub2ind(B, sequences(v,ix), sequences(b,ix), ix))));
        end
        children{a} = find(tree_(:,1) == a);
        children{b} = [children{b}; v];
        ev_b = evidence(b);  % ev_b = true iff b was not an evidence node before
        evidence(b) = true;
        evidence(a) = max([evidence(children{a}); ~isempty(find(t==a, 1))]);     
        bias = length(pool)/(length(pool) - ~evidence(a) + ~ev_b);        
    end
    [~,rev_transition_prob] = manipulate(tree_, code, children);

    acceptance_prob = bias*exp(ll_diff + log(rev_transition_prob) - log(transition_prob));
    accepted = rand<acceptance_prob;
    if accepted || (greedy == -1 && ll_diff > 0)
        codes(code(3)) = codes(code(3))+1;
        tree = tree_;
        total_ll_diff = total_ll_diff + ll_diff;
    else
        if a ~= b
            children{b} = children{b}(1:end-1);
            children{a} = [children{a};  v];
            evidence(a) = true;
            evidence(b) = ev_b;
        end
        codes(5) = codes(5)+1;
    end    
end
assert(isequal(evidence,nodes_with_evidence(tree,t)));

end


function [tree transition_prob code] = manipulate(tree, code, children) 

unstable = false;

% v is given in code(2).  
assert(code(1) ~= 0);
v = code(1);

assert(tree(1) == 0); % first node is the root

% let a be its parent cell.
a = tree(v,1);
assert(a~=0); % don't handle root

% let t2 be the event where v was born. 
t2 = tree(v,2);

% let t1 and t3 be the previous and next events associated with a. (always
% exist)

[t3 t1 next prev] = adjacent_events(tree,v, children);

% let t0 be the previous event going up the tree from t1, 
% let t0_ be the next event going down the *alternative* path from t1.
[t0_ t0 prev_next] = adjacent_events(tree,prev, children);
if prev_next == v
    t0_ = first_event(tree,prev,children);   
end

% let t3_ be the next event associated with t2.
% the new time of event t2 must be before event t3_.
t3_ = first_event(tree,v, children);

if t3_ < t3
    t3 = t3_;
    t4 = -1; t4_ = -1;
else
% let t4 be the next event associated with a going down the tree from t3 
% let t4_ be the next event associated with the new cell born in t3
    t4 = adjacent_events(tree, next, children);
    t4_ = first_event(tree, next, children);
end

% some of these other events may not exist.
% choose one of the existing events, uniformly.
z = double([t0 t0_ t4 t4_] >= 0);  % z is binary. z==0 if that option is disabled.
if all(z == 0)
        z(2) = 1;
        t0_ = t1;
end
[k transition_prob] = lmnrnd(log(z),code(2));
if z(k)<= 0
    fprintf('Something weird is going on\n');
    code
    z
    k
    keyboard;
    % What happens here is that there are two cells with exactly the same
    % birth time, and the current node is "between them".  But because the
    % times are the same, it did not change the ordering of the events so
    % some things are messed up.  DONE: make sure no events are too close.
end

if ~unstable && code(2) ~= 0, return; end % reverse move - no need for more computation
        

% assume t0 was chosen (same idea for t0_).
% select a new place/time for the event t2 by sliding it on the track 
%     t0---t1---t2---t3
if k == 1 % backward-backward
%    z = [t1 min(t3,t3_); t0 t1];
    z = [t1 t3; t0 t1];
    pr = [a tree(prev,1)];
    k_rev = 3+(pr(2) ~= a);
elseif k==2 % 
    z = [t1 t3; t1 min(t0_,t3_)];
    if prev == a
        pr = [a tree(prev,1)];
    else
        pr = [a prev];
    end
    k_rev = 2;
elseif k==3 % forward-forward1
% assume t4 was chosen (same idea for t0_).
% select a new place/time for the event t2 by sliding it on the track 
%     t1---t2---t3---t4
% again, the new time of event t2 must be before event t3_.
    z = [t1 t3; t3 min(t3_,t4)];
    pr = [a a];   
    k_rev = 1;
else % k==4
    z = [t1 t3; t3 min(t3_,t4_)];
    pr = [a next];
    k_rev = 1;
end


% new addition - don't allow to stay in the current interval!
y = z(:,2)-z(:,1);
if unstable && y(2) > 0
    y(1) = 0;
    transition_prob = transition_prob/y(2);  % unstable - move out of own region
else 
    assert(y(1) > 0);
end
if code(2) ~= 0, return; end  % reverse move - no need for more computation

t = sum(y)*rand;
assert(t>=0);

l = 1;
if t>y(1)
    l = 2;
    t = t-y(1);
end
t = z(l,1)+t;

%    v was born at time t and its parent is pr(l) 
tree(v,1) = pr(l);
tree(v,2) = t;

if tree(v,2)>=tree(v,3), fprintf('Something wierd is happening\n'); keyboard; end

% To do the reverse move - choose the same cell v.
% choose t0' = t1, and the track is identical to 
%     t0'--t1'--t2'--t3'
%     t1---t2---t3---t4
% 
if l==1, k_rev = k; end
code = [v k_rev k];

% fprintf('sum(y) = %.2f  k=%d k_rev=%d, z=[%.2f %.2f; %.2f %.2f]  t2=%.2f t2_new=%.2f ', sum(y), k, k_rev, z(1,1), z(1,2), z(2,1), z(2,2), t2, t);
assert(k ~= 0);
assert(v ~= 0);

end


function [t_next t_prev next prev ] = adjacent_events(tree,v, children)
    prev = 0; next = 0; t_prev=-1; t_next=-1;
    if v <= 0, return; end
    a = tree(v,1);
    if a == 0, return; end
    ix = children{a};
%     ix_ = find(tree(:,1) == a);
%     assert(isequal(sort(ix(:)),sort(ix_(:))));
    
    % sort the births according to time.
    [~,o] = sort(tree(ix,2));
    ix = ix(o);
    
    i = find(ix == v);
    if i == 1
        prev = a;
    else
        prev = ix(i-1);
    end
    t_prev = tree(prev,2);
    
    if i == length(ix)
        next = 0; %-a;
        t_next = tree(a,3);  % when a is dead
    else
        next = ix(i+1);
        t_next = tree(next,2);
    end  
    

end

function [t first] = first_event(tree,v, children)
    t = -1; first = 0;
    if v == 0, return; end
    ix = children{v};
%     ix_ = find(tree(:,1) == v);
%     assert(isequal(sort(ix(:)),sort(ix_(:))));
    
    if isempty(ix)
        first = 0;
        t = tree(v,3);
    else 
        [~,i] = min(tree(ix,2));
        first = ix(i);
        t = tree(first,2);
    end
end


function sanity(tree)
    for u=2:size(tree,1);
        assert(tree(u,2)>tree(tree(u,1),2));
        assert(tree(u,2)<=tree(tree(u,1),3));
        assert(tree(u,2)<=tree(u,3));
    end
end




