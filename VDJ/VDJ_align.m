function [L_aligned R R_aligned trimmed_VJ st al] = VDJ_align(reads, V_seq, J_seq, edge_VJ)
% edge_VJ is how much from the V and J genes you want to use inside R

%  Trim the edges of the reads so they start and end in perfect alignment
cut_v = length(V_seq)-edge_VJ(1)+1;
cut_j = length(J_seq)-edge_VJ(2)+1;
[L_aligned R st1] = trim_and_fix_reads(reads, V_seq, 0, cut_v);
if any(st1(2,:) > 6), 
    fprintf('%d letters were erased from beginning of some reads!\n', max(st1(2,:)));
end
if any(st1(1,:) < st1(2,:))
    fprintf('Some reads end %d letters before the J starts!\n', min(st1(1,:) - st1(2,:)));    
end    

%[L_aligned R] = trim_and_fix_reads(R, V_seq(cut_v:end), 0, 20);    
[R_aligned R st2] = trim_and_fix_reads(R, J_seq, 1, cut_j);
if any(st2(2,:) > 16), 
    fprintf('%d letters were erased from ending of some reads!\n', max(st2(2,:)));
end
if any(st2(1,:) < st2(2,:))
    fprintf('Some reads end %d letters after the J ends!\n', min(st2(1,:) - st2(2,:)));    
end    

st = [st1; st2];

trimmed_VJ = [0 0];
R_aligned = int16(cell2mat(R_aligned));
L_aligned = int16(cell2mat(L_aligned));

trimmed_VJ(1) = find(diff(sum(L_aligned == 5,1))/sum(diff(sum(L_aligned == 5,1)))>0.1,1, 'last');
L_aligned = L_aligned(:,trimmed_VJ(1)+1:end);

trmJ = find(diff(sum(R_aligned == 5,1))/sum(diff(sum(R_aligned == 5,1)))>0.1,1, 'first');
if ~isempty(trmJ)
    trimmed_VJ(2) = size(R_aligned,2) - trmJ;
end
R_aligned = R_aligned(:,1:end-trimmed_VJ(2));


[~, d, al] = profileHMM_align(R, {V_seq(end-edge_V+1:end)}, {rep.D.Sequence}, {J_seq(1:edge_J)}, noise);



end