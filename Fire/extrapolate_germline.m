function germline = extrapolate_germline(consensus, edge_VJ, trimmed_VJ, rep)

noise = 0.001;
map = []; map('ACGTNR') = [1:5 5];
dict = 'ACGTN';

V_seq = rep.V.Sequence;
J_seq = rep.J.Sequence;

[~, d, al] = profileHMM_align({consensus}, {V_seq(end-edge_VJ(1)+1:end)}, {rep.D.Sequence}, {J_seq(1:edge_VJ(2))}, noise);
G = map([V_seq(trimmed_VJ(1)+1:end-edge_VJ(1)) dict(al.germline) J_seq(edge_VJ(2)+1:end-trimmed_VJ(2))]);

lastV  = length(V_seq)-trimmed_VJ(1)-al.eaten(1);
firstJ = length(J_seq)-trimmed_VJ(2)-al.eaten(4);
firstJ = length(G)+1 - firstJ; % correct to index from left
germline.in1 = lastV  + ([1:length(al.N1)]);
germline.in2 = firstJ - ([length(al.N2):-1:1]); 
eaten = [trimmed_VJ(1) al.eaten trimmed_VJ(2)];

% construct germline header V_J_D_eaten_lastV_firstJ_N1_N2
germline.header = ['germline_' rep.V.Header '_' rep.J.Header];
if d ~= 0 % D is known
    germline.header = [germline.header '_' rep.D(d).Header sprintf('_%d', eaten)];
else
    germline.header = [germline.header '_IGHD-Unknown' sprintf('_%d', eaten)];
end
germline.header = [germline.header sprintf('_%d', [lastV firstJ])];

germline.seq = G;
germline.Sequence = dict(G);    
germline.eaten = eaten;
germline.d = d;
germline.lastV_firstJ = [lastV firstJ];
germline.frame_shift = mod(eaten(1)+eaten(6) + length(germline.seq) - 1, 3);
germline.n1 = dict(al.N1);
germline.n2 = dict(al.N2);

