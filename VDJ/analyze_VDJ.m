%function [v,j,err, d, al] = analyze_VDJ(seq, rep, v, j, edge_VJ)
% given the sequence seq annotates it with V,D,J combination. 
% v,d,j are the indices in the repertoire referencing the inferred v,d,j.
% V and J are classified independently, and then we go over all possible
% Ds and infer the deletions and additions in the junctions.
function [v,j,err, d, al] = analyze_VDJ(seq, rep, v, j, edge_VJ)

if length(seq)-1 <= 150
    v=-1; d=-1; j=-1; al=''; err= 10;
    fprintf('analyze_VDJ: input sequence is too short (%d)\n', length(seq))
    return;
end
d = 0; 
if nargin < 3, v = 0; end
if nargin < 4, j = 0; end
if nargin < 5, edge_VJ = [40 25]; end

V = []; N1 = []; D = []; N2 = []; J = [];
% Classify V

vst = ones(length(rep.V), 2);
if v == 0
    sc = zeros(length(rep.V), 1);
    for v = 1:length(rep.V)
            [~,~,vst(v,:),sc(v)] = trim_and_fix_reads({seq}, rep.V(v).Sequence, 0);        
    end
    [~,v] = max(sc);    
end
% note taht trim_and_fix_reads passes V in before seq to nwalign, so
% thus vst is from swalign(V, seq)
% vst(v,:)=[dropped prefix from rep, dropped prefix from seq]
J_start_search = 1 + vst(v,2) - vst(v,1) + length(rep.V(v).Sequence) - edge_VJ(1);

J_search_length = length(seq) - J_start_search; % the subseq in which we search for J -- how many bases does it contain?
if J_search_length < 10 || J_start_search < 50
    d=-1; j=-1; al=''; err= 11;
    fprintf('analyze_VDJ: either J search length (%d) or start search point of J (%d) is too early \n', J_search_length, J_start_search)
    return
end

% Classify J
if j == 0
    sc = zeros(length(rep.J), 1);
    for j=1:length(rep.J)
            [~,~,~,sc(j)] = trim_and_fix_reads({seq(J_start_search:end)}, rep.J(j).Sequence, 1);            
    end
    [~,j] = max(sc);
end

if nargout <=2, return; end

V_seq = rep.V(v).Sequence;
J_seq = rep.J(j).Sequence;
[L_aligned R R_aligned trimmed_VJ st err] = VJ_align({seq}, V_seq, J_seq, edge_VJ);

if nargout <=3, return; end

% Classify D and annotate
[~, d, al] = profileHMM_align(R, ...
             {V_seq(end-edge_VJ(1)+1:end)}, {rep.D.Sequence}, {J_seq(1:edge_VJ(2))}, 0.05);


% output         
map('ACGTN') = int8(1:5);
V_seq = map(V_seq);
J_seq = map(J_seq);

al.V  = int8([L_aligned al.V]);
al.V_ = int8([V_seq(trimmed_VJ(1)+1 : end-edge_VJ(1)) al.V_]);

al.J  = [al.J R_aligned ];
al.J_ = [al.J_ J_seq(edge_VJ(2)+1 : end-trimmed_VJ(2)) ];

al.germline = [V_seq(trimmed_VJ(1)+1 : end-edge_VJ(1)) al.germline J_seq(edge_VJ(2)+1 : end-trimmed_VJ(2))];

al.N1 = int8(al.N1);
al.N2 = int8(al.N2);
al.D  = int8(al.D);
al.D_ = int8(al.D_);
al.trimmed_VJ = trimmed_VJ;

end
