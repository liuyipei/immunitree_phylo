function [dg dr ] = analyze_intra_read_similarity(clone_file, inter_replica, is_variant)    

    dg = zeros(5,2);
    dr = zeros(5,2);

    [clone germline id src replicas] =  play_with_clone(clone_file);
    if is_variant
        clone = [germline.seq; mode(clone,1)];  % returns germ_to_consensus instead of mean germ_to_read
    else
        clone = [germline.seq; clone];
    end
    
    % calcualte D, distance matrix between reads. Exclude mismatches with gaps.
    B = pdist(double(clone), 'hamming')*size(clone,2);
    R = pdist(double(clone == 5), 'hamming')*size(clone,2);
    D = B-R;
    if nargout == 0, imagesc(D); return; end
    
    % calculate mean and std of distance to germline
    dist_to_germline = D(1:size(clone,1)-1);
    dg(5,1) = mean(dist_to_germline);
    dg(5,2) = std(dist_to_germline);

    % calculate mean and std of pairwise read distance
    dist_intra_read = D(size(clone,1):end);
    D = squareform(dist_intra_read, 'tomatrix');
    N = size(D,1);

    % if we want to remove inter-replica bias
    if inter_replica
        % put NaN in all intra-replica pairs.
%        replicas = [data.subject_num(id) data.FR(id) data.rep(id)-'a'+1];
        [~,~,tag] = unique(replicas, 'rows');
        for s=1:max(tag)
            ix = tag == s;
            D(ix,ix) = NaN;            
        end
        D(sub2ind([N N], 1:N, 1:N)) = 0 ; % set diagonal to 0
        dist_intra_read = squareform(D, 'tovector');    
        dist_intra_read = dist_intra_read(~isnan(dist_intra_read));
    end
    
    dr(5,1) = mean(dist_intra_read); 
    dr(5,2) = var(dist_intra_read); 
    
    
    % Dsrc = zeros(4);
    
    % (same per source)
    if ~isempty(D)
        for s = 1:4
            ix = src == s;

            % distance to germline of source
            dg(s,1) = mean(dist_to_germline(ix));
            dg(s,2) = std(dist_to_germline(ix));

            tmp = squareform(D(ix,ix), 'tovector');
            tmp = tmp(~isnan(tmp));
            dr(s,1) = mean(tmp);
            dr(s,2) = std(tmp); 
    %         for s_ = s+1:4
    %             iy = src == s_;
    %             Dsrc(s_,s) = mean(reshape(D(ix,iy), [],1));
    %         end
        end 
    %    Dsrc = squareform(Dsrc);
    end
    
end


function test()
%%
[dg dr Dsrc] = analyze_intra_read_similarity([tree_dir 'all_clones/' trees(9).name])
end