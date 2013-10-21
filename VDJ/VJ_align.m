% alignes reads to V and J (not D) and trims the edges if needed.
function [L_aligned R R_aligned trimmed_VJ st err] = VJ_align(reads, V_seq, J_seq, edge_VJ)
% edge_VJ is how much from the V and J genes you want to use inside R
err = false;

%  Trim the edges of the reads so they start and end in perfect alignment
cut_v = length(V_seq)-edge_VJ(1)+1;
cut_j = length(J_seq)-edge_VJ(2)+1;
[L_aligned R st1] = trim_and_fix_reads(reads, V_seq, 0, cut_v);
if any(st1(2,:) > 2) % remember that the reads were passed in a a set
    % err = 1;
    % fprintf('VJ_Align: Some prefix of up to %d letters of the read were unaligned to any V.\n', max(st1(1,:)));
end
if any(st1(1,:) > 2) % remember that the reads were passed in a a set
    % err = 20;
    % fprintf('VJ_Align: Some prefix of up to %d letters of the repertoire V were unaligned.\n', max(st1(1,:)));
end

R_lengths = cellfun(@length, R);
if any(R_lengths<30) % remeber that the reads were passed in a a set
    err = 21;
    fprintf('VJ_Align: Some sequence following aligned V is too short.\n');
    fprintf(',%d', R_lengths);
    R_aligned = int16([]); st = [0,0];
    trimmed_VJ =[0 0];
    return
end

[R_aligned R st2] = trim_and_fix_reads(R, J_seq, 1, cut_j);
if any(st2(2,:) > 5)%2)
    % fprintf('VJ_Align Warn: Some suffix of up to %d letters of the read were unaligned to any J.\n', max(st1(1,:)));
    %err = 3;
end
if any(st2(1,:) > 5)%2)
    % err = 25;
    % fprintf('VJ_Align Warn: Some suffix of up to %d letters of the read were unaligned to its best rep J.\n', max(st2(1,:)));
end

if any(cellfun(@length, R) == 0)
    err = 5;
    fprintf('VJ_Align Err: No variable region left!\n');
end
st = [st1; st2];

trimmed_VJ = [0 0];
R_aligned = int16(cell2mat(R_aligned)); 
L_aligned = int16(cell2mat(L_aligned));

orig_joni_trim_qc = false;
if orig_joni_trim_qc
    % trim V one base before the first location 90% of the reads are supported
    non_gaps = sum(L_aligned ~= 5,1) / size(reads,1);
    trmV = find(non_gaps > 0.9, 1, 'first');
    if isempty(trmV)
        fprintf('VJ_align Err: Garbage read, poor V alignment\n'); 
        err = 8;
        return;
    end
    trimmed_VJ(1) = trmV-1;
    L_aligned = L_aligned(:,trimmed_VJ(1)+1:end);
    trmV = find(non_gaps > 0.9, 1, 'last'); % last location where there is no gap
    if size(L_aligned,2)-trmV > 0
        fprintf('VJ_align Err: Too much of the V on the junction side was eaten.\n');
        err = 6;
    end

    % trim one base after the *last* location 90% of the reads are supported
    non_gaps = sum(R_aligned ~= 5,1) / size(reads,1);
    trmJ = find(non_gaps > 0.9, 1, 'last');
    if isempty(trmJ)
        fprintf('VJ_align Err: garbage read, poor J alignment\n'); 
        err = 9;
        return;
    end
    trimmed_VJ(2) = size(R_aligned,2) - trmJ;
    R_aligned = R_aligned(:,1:end-trimmed_VJ(2));
    trmJ = find(non_gaps > 0.9, 1, 'first'); % first location where there is no gap
    if trmJ > 1
        fprintf('VJ_align Err: Too much of the J on the junction side was eaten.\n');
        err = 7;
    end
end
end
