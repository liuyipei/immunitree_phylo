function triplets_based_statistics
%%
clone_triplets =  struct('chainy', -1, 'bushy', -1);
clone_triplets = clone_triplets(ones(1,N), 1);

for i=1:N
    if mod(i,1000)==0, fprintf('%d ', i); end
    [clone germline id src Z] =  play_with_clone([tree_dir trees(i).name]);
    
    % prepare list of triplets
    if size(clone,1) == 3
        triplets = [1 2 3];
    elseif size(clone,1) == 4 
        triplets = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
    elseif size(clone,1) == 5
        triplets = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
        triplets = [triplets ; 1+triplets(2:end,:) ; 1 4 5; 1 2 5; 1 3 5];
    else
        triplets = ceil(size(clone,1)*rand(100,1));
        triplets = [triplets ceil((size(clone,1)-1)*rand(100,1))];
        triplets = [triplets ceil((size(clone,1)-2)*rand(100,1))];
        ix = triplets(:,2) >= triplets(:,1);
        triplets(ix,2) = triplets(ix,2)+1;
        ix = triplets(:,3) >= triplets(:,2);
        triplets(ix,3) = triplets(ix,3)+1;        
    end
    
    % calcualte D, distance matrix between reads. Exclude mismatches with gaps.
    B = pdist(double(clone), 'hamming')*size(clone,2);
    R = pdist(double(clone == 5), 'hamming')*size(clone,2);
    D = B-0.5*R;    % the 0.5 is a necessary hack to make sure that the 
                    % triangle inequility still holds
    D = squareform(D);
    
    % given 3 reads, x, y, z
    % For each read, find if it is in the path between the other two by
    % computing for each read y U(y) = D(x,y)+D(y,z)-D(x,z).    
    V = [1 1 -1; 1 -1 1; -1 1 1];
    ix =     sub2ind(size(D), triplets(:,1), triplets(:,2));
    ix = [ix sub2ind(size(D), triplets(:,2), triplets(:,3))];
    ix = [ix sub2ind(size(D), triplets(:,1), triplets(:,3))];
    U = D(ix)*V;    
    
    % the bushiness of a triplet x y z is min(U(x) U(y) U(z))
    % the bushiness of a clone is the average bushiness of triplets
    clone_triplets(i).bushy = 0.5*mean(min(U,[],2));
    
    % another approach:  Find how often one read is between the other two
    % on the tree.
    
    % Construct D.  D(i,j) = true <==> j is ancestor of i,
    T = Z.a.tree;
    M = size(T,1);
    [anc, o, d] = order_tree(Z.a.tree);
    D = false(M,M);
    for j=1:length(o)
        D(o(j), anc(j,1:d(j))) = true;
    end
    
    
    cnt = 0;
    for j=1:size(triplets,1)
        % trp indicates the tree nodes of the reads in the current triplet
        trp = false(1,M);    
        trp(Z.a.t(triplets(j,:))) = true;
        
        % if all the nodes in trp appear as ancestors anywhere ==> cnt++
        ix = sum(trp(ones(M,1),:) <= D,2);
        cnt = cnt + any(ix == M);
    end
%    cnt/size(triplets,1)
    
    clone_triplets(i).chainy = cnt/size(triplets,1);
    
end
fprintf('\n');

%%

labels = label_reads('source', tree_dir);
X = [ [clone_triplets.chainy]' [clone_triplets.bushy]'];
%X(:,1) = X(:,1)/100;

plot_multi_hist_with_significance(X, labels, {'chainy', 'bushy'}, [0.1 1]);
saveas(gcf, '~/www/research/donors/bushy_chainy.jpg');
saveas(gcf, '~/www/research/donors/bushy_chainy.fig');


end

