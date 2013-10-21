% how confidence are we in the tree results?

K = 10;
J = 10;
stats = zeros(K,J,7);
%%
for k=1:K
    %% create synthetic example
        me = 0.5; 
        alpha = 0.5;
        sigma = 0.025;
        epsilon = 0.003;
        N = 2000;
        L = 400;
        [reads, sequences, T, t] = ...
                generate_data(N, L, me, alpha, epsilon, sigma);        
        view_tree(T);
    %%
    %save('synthetic tree.mat', 'reads', 'sequences', 'T', 't', 'me', 'alpha', 'sigma', 'epsilon', 'N', 'L');

    %% run J times
    for j=1:J
        params = [];
        params.me = 0.5;
        params.alpha = 0.5;
        params.epsilon = 0.001;
        params.sigma = 0.025;
        params.greedy = 0;
        [T_ seqs_ t_]= DPtrees_inference_binary(reads, 50, params);

    %% Turn into a graph
        UG = get_graph_from_tree(T, true, true);
        D = graphallshortestpaths(UG, 'directed', false)

        UG_ = get_graph_from_tree(T_, true, true);
        D_ = graphallshortestpaths(UG_, 'directed', false)

        B = D(t,t);
        B_ = D_(t_,t_);
        % Take a snapshot after every X iterations

    %% consider the average distance between two reads.
        figure(1); rand_index = confusion_matrix(t, t_)
        max_dist = max(max(D(:)), max(max(D_(:))));
        B = squareform(D(t,t));
        B_ = squareform(D_(t_,t_));
        Z = B_-B;

        total_pairs = N*(N-1)/2;
        fprintf('There is a total of %d total pairs.\n', total_pairs);

        were_together = sum(B==0);
        fprintf('Of these pairs, %.1f%% are in the same class in truth.\n', ...
            100*were_together/total_pairs);
        
        % how many pairs of reads where in distance 0 before?  Where are they now?
        still_together = sum(B == 0 & B_ == 0);
        fprintf('Of these pairs, %.1f%% are still together now, which is %.1f%% of the total pairs.\n',...
            100*still_together/were_together, 100*still_together/total_pairs);

        % no. of times dist hasn't changed.
        same_dist = sum(Z == 0);
        fprintf('%.1f%% of all pairs have same tree distance before and after.\n',...
            100*same_dist/total_pairs);
        
        avg_dist_change = sum(abs(Z))/(total_pairs - same_dist);
        tendency = sum(Z)/(total_pairs - same_dist);
        fprintf('Those pairs that changed their distance, changed it by %.1f edges on average.\n',...
            avg_dist_change);
        fprintf('Overall, the distances are now %.1f edges further.\n', tendency);
        figure(2); hist(Z, -max_dist:max_dist);
        stats(k,j,:) = [total_pairs were_together still_together same_dist avg_dist_change tendency rand_index];
    end

end



%%

% how many reads where in distance 1 before?  Where are they now?

% Look a