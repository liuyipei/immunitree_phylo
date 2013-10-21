function st = label_reads(label_str, tree_dir, param)

if nargin < 2
    tree_dir = '/afs/cs/u/joni/scratch/data/lymph/t5.fa_results/';
end

load([tree_dir 'all_clones.mat'], 'clone_ids', 'clone_src_');

names = [];

if ~isstr(label_str) % custom labels
    labels = label_str;
    label_str = 'label'; 
    max_label = max(labels);
    
elseif strcmp(label_str, 'source')
    [~,labels]  = max(clone_src_,[], 2);

    ix = sum(clone_src_ > 0, 2) > 1;
    labels(ix) = 5;    
    max_label = 5;
    names = {'blood', 'spleen', 'mes', 'med', 'multi'};
    
elseif strcmp(label_str, 'pairwise only')
    labels = zeros(size(clone_src_,1),1);
    
    cur_label = 1;
    for i=1:4
        for j=(i+1):4
            labels(clone_src_(:,i)>1 & clone_src_(:,j)>1) = cur_label;
            cur_label = cur_label+1;
        end
    end
    ix = sum(clone_src_ > 1, 2) > 2;
    labels(ix) = cur_label;
    max_label = cur_label;
    names = {'12', '13', '14', '23', '24', '34', 'triplets'};
    
    % at this point all the clones with exactly two sources with > 1 read 
    % are label.  Clones with 3 sources of more with >1 reads are labeled 
    % as well.  Other clones (only one source with >1 reads) are 0.
    
elseif strcmp(label_str, 'pairwise source')
    [~,labels]  = max(clone_src_,[], 2);

    ix = sum(clone_src_ > 0, 2) > 1;
    labels(ix) = 0;
    
    % at this point all the pure clones are labeled, and all the non-pure
    % clones have label 0.
    
    cur_label = 5;
    for i=1:4
        for j=(i+1):4
            labels(clone_src_(:,i)>1 & clone_src_(:,j)>1) = cur_label;
            cur_label = cur_label+1;
        end
    end
    ix = sum(clone_src_ > 1, 2) > 2;
    labels(ix) = cur_label;
    max_label = cur_label;
    names = {'1', '2', '3', '4', '12', '13', '14', '23', '24', '34', 'triplets'};
    
    % at this point all the pure clones are labeled and all the clones with 
    % exactly two sources with > 1 read are label.  Clones with 3 sources
    % of more with >1 reads are labeled as well.  Other clones (not pure
    % but with only one source with >1 reads) are 0.
    
elseif strcmp(label_str, 'source organ')
    [~,labels]  = max(clone_src_,[], 2);
    labels(labels == 4) = 3;
    
    ix = sum(clone_src_ > 0, 2) > 1;
    labels(ix) = 4;
    max_label = 4;
    names = {'blood', 'spleen', 'LNs', 'multi'};
    
else
    if ~strcmp(label_str, 'default')
        fprintf('unrecognized label, using default\n');
    end
    labels = ones(N,1);
    max_label = 1;
end

st.l = labels;
st.num = max_label;
if isempty(names), names = strread(sprintf('%d\t',1:max_label), '%s'); end
st.names = names;
st.str = label_str;
st.patients = clone_ids(:,1);

end