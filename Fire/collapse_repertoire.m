function [rep newrep data] = collapse_repertoire(rep, keep, data)

if nargin<2 keep = 1:length(rep.V); end

rep_ = rep;

rep.V = rep.V(keep);
genes = {rep.V.Header};

% add '*' to the end of the genes to that genes like IGHV4-VH4*80P are 
% included in this collapse.  If it is commented, then these genes will 
% not be included in the new repertoire
% genes = cellfun(@(x) [x '*'], genes, 'uniformoutput', false);

% trim for each gene what comes before the first '*'
genes = cellfun(@(x) x(1:find(x == '*',1)-1), genes, 'uniformoutput', false);

% show consolidated genes
[genes I V] = unique(genes);
for i=1:length(genes), fprintf('%d: %s\n',i, genes{i}); end


V_ = (length(I)+1)*ones(1,length(keep));
V_(keep) = V;
V = V_;

rep_.col_V = genes;
rep_.map_V = V;

% do the same for J

genes = {rep.J.Header};
genes = cellfun(@(x) [x '*'], genes, 'uniformoutput', false);
genes = cellfun(@(x) x(1:find(x == '*',1)-1), genes, 'uniformoutput', false);
[genes I J] = unique(genes);
for i=1:length(genes), fprintf('%d: %s\n',i, genes{i}); end

rep_.col_J = genes;
rep_.map_J = J;

rep = rep_;


if nargout > 1
    assert(exist('data', 'var') == 1);
    newcode = [rep.code '_collapsed'];

    %  Create a new reportoire by taking the most commong gene in each group
    newrep = [];
    [~, ~, big_hist] = collect_VJs(data, get_reads(data,1:40), rep.code, 0);
    sumV = sum(big_hist,2)';
    sumJ = sum(big_hist);

    % V Genes
    fprintf('\n');
    lens = cellfun(@length, {rep.V.Sequence});      
    %TODO: do not include Vs that have only 1-2 reads aligned to them
    % IGHV5-a*02 - I'm looking at you!
    for i=1:length(rep.col_V)    
        ix = find(rep.map_V == i);
        if max(sumV(ix))>2, ix = ix(sumV(ix)>2); end        
        [~, o] = sortrows([lens(ix)' sumV(ix)'], [-1 -2]);
        j = ix(o(1));
        newrep.V(i).Header = rep.col_V{i};
        newrep.V(i).Sequence = rep.V(j).Sequence;
        fprintf('%2d: %s \t--> %s \t[alleles: %2d;  reads: %d/%d]\n', ...  % 'len: %d/%d)%s\n', 
            i, rep.col_V{i}, rep.V(j).Header(1:min(end,16)), length(ix), ...
            sumV(j), sum(sumV(ix))); %, lens(j), max(lens(ix)), as);
    end

    newrep.V = newrep.V';
    newrep.D = rep.D;

    % J Genes
    fprintf('\n');
    lens = cellfun(@length, {rep.J.Sequence});      
    for i=1:length(rep.col_J)    
        ix = find(rep.map_J == i);
        [~, o] = sortrows([lens(ix)' sumJ(ix)'], [-1 -2]);
        j = ix(o(1));
        newrep.J(i).Header = rep.col_J{i};
        newrep.J(i).Sequence = rep.J(j).Sequence;
        fprintf('%2d: %s \t--> %s \t[alleles: %2d;  reads: %d/%d]\n', i, rep.col_J{i}, ...
            rep.J(j).Header(1:min(end,16)), length(ix), sumJ(j), sum(sumJ(ix)));
    end
    newrep.J = newrep.J';
    newrep.code = newcode;

end

if nargout > 2 % re-map the reads to the collapsed repertoire
    data.(newcode) = data.(rep.code);
    
    % V genes
    iv = [1 4:size(data.(rep.code),2)];
    ihmmune = data.(rep.code)(:,iv);
    ix = ihmmune > 0;
    ihmmune(ix) = rep.map_V(ihmmune(ix));
    ix = ihmmune == length(newrep.V)+1;
    fprintf('%d reads mapped to pseudo genes, these mappings are removed.\n', sum(ix(:)));
    ihmmune(ix) = 0;
    data.(newcode)(:,iv) = ihmmune;
        
    % J genes
    ix = data.(rep.code)(:,3) > 0;
    data.(newcode)(ix,3) = rep.map_J(data.(newcode)(ix,3));
end

end



function show_distance_between_repoertoire_genes(dist_bu)
%%
    figure; 
    [sorted_J o] = sort(J);
    dist = dist_bu;
    dist_ = squareform(dist);
    %dist_ = dist_(keep,keep);
    dist_ = dist_(o,o);
    imagesc(dist_, [0 1]);
    colormap('gray');
    plot_class_lines(sorted_J, true, 'm');
    plot_class_lines(sorted_J, false, 'm');
%%
    figure;
%    [sorted_V oV] = sort(rep.map_V);
    lens = cellfun(@length, {rep.V(oV).Sequence});      
    bar(lens);
    plot_class_lines(sorted_V, false, 'm');
    ylim([150 350]);
    title('length of V genes in repertoire')    
end
