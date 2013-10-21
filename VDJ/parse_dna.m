function [Z stamps germline] = parse_dna(samples, FRs)
%% 

% Chose which sample:
ix_stamps = 0;
Z = [];
stamps = [];
base_dir = '/afs/cs/u/joni/scratch/data/Fire_FL/';
FR_dir = {'S57to73FR1', 'S57to73FR2'};
if nargin<2, FRs = 1:length(FR_dir); end
FRs
M = length(FRs);
for nSample = samples
    Z_ = cell(M, 1);
    stamps_ = Z_;
    for k=1:M
        files = dir(sprintf('%s/%s/S%d*', base_dir, FR_dir{FRs(k)}, nSample));
        assert(length(files) == 1);
        filename = [base_dir FR_dir{FRs(k)} '/' files.name]
        [Z_{k} stamps_{k}]= parse_dna_file(filename, true);
        stamps_{k} = stamps_{k}+ix_stamps;
        ix_stamps = max(stamps_{k});
    end
    trunc = min(cellfun(@(x) size(x,2), Z_));
    Z = [Z; cell2mat(cellfun(@(x) x(:,1:trunc), Z_, 'UniformOutput', false))];
    stamps = [stamps cell2mat(stamps_')];
end

clear Z_ stamps_

cons_dist = zeros(ix_stamps);
for i=1:ix_stamps
 for j=i+1:ix_stamps
	seq1 = Z(find(stamps==i,1), :);
	seq2 = Z(find(stamps==j,1), :);
	cons_dist(i,j) = sum(seq1 ~= seq2);
 end
end
cons_dist

%  Load repertoire of D's and the chosen V and J genes
    rep = load_repertoire('ihmmune');
%     rep.V = fastaread('V.fa');
%     rep.J = fastaread('J.fa');
%     rep.D = fastaread('D.fa');

%  Prepare data
%Z = seqrcomplement(Z);
Z = fliplr(Z);
Z = 5-Z;
Z(Z == 0) = 5;
cons = Z(1,:);
%cons = seqrcomplement(cons);

% for each V-J
sc = zeros(length(rep.V), 1);
for v = 1:length(rep.V)
        [~,~,~,sc(v)] = trim_and_fix_reads({cons}, rep.V(v).Sequence, 0);        
end
[~,v] = max(sc);
[~,~,st] = trim_and_fix_reads({cons}, rep.V(v).Sequence, 0);
if (st(2) > 6), assert(false); end
V_seq = rep.V(v).Sequence((st(1)-st(2)+1):end);
%V_seq = rep.V(v).Sequence(st(1):end);
lenV = length(V_seq);
germline.frame_shift = mod(st(1)-st(2), 3);
%germline.frame_shift = mod(st(1)-1, 3);

sc = zeros(length(rep.J), 1);
for j=1:length(rep.J)
        [~,~,~,sc(j)] = trim_and_fix_reads({cons(lenV-40:end)}, rep.J(j).Sequence, 1);            
end
[~,j] = max(sc);
[~,~,st] = trim_and_fix_reads({cons(lenV-40:end)}, rep.J(j).Sequence, 1);
if (st(2) > 6), assert(false); end
J_seq = rep.J(j).Sequence(1:(end-st(1)+st(2)) );
% assert(st(2) == 1);
% J_seq = rep.J(j).Sequence(1:(end-st(1)+1) );


%  Map reads to profile HMM
fprintf('Mapping the reads to profile HMMs...\n');
% tic;
% [profileHMM_score profileHMM_id] = profileHMM_align({cons}, {V_seq}, {rep.D.Sequence}, {J_seq});
% toc;    

% align to profile HMMs
sc = zeros(length(rep.D), 1);
params = profileHMM_get_params();
params.sigma = 0.05;
for d = 1:length(rep.D)
    [~,sc(d)] = ...
        profileHMM_annotate(cons, V_seq, rep.D(d).Sequence, J_seq, params);
end
[~,d] = max(sc);
a = profileHMM_annotate(cons, V_seq, rep.D(d).Sequence, J_seq, params);
[V,N1,D,N2,J, eaten] = deal(a.V, a.N1, a.D, a.N2, a.J, a.eaten);
germline.al = a;
germline.vdj = [v d j];

map('ACGTNacgtn') = [1:5 1:5];
D_seq = map(rep.D(d).Sequence((eaten(2)+1):(end-eaten(3))));
V_seq = map(V_seq(1:length(V)));
J_seq = map(J_seq(end-length(J)+1:end));
germline.seq = [V_seq 5*ones(size(N1)) D_seq 5*ones(size(N2)) J_seq];
dict = 'ACGT-';
fprintf('%s\n', dict(germline.seq));
fprintf('%s\n%s\n%s\n', rep.V(v).Header, rep.D(d).Header, rep.J(j).Header);

dist_to_cons = sum(cons(ones(1,size(Z,1)),:) ~= Z, 2)-sum(Z==5, 2)-sum(cons==5, 2);
dist_to_germ = sum(germline.seq(ones(1,size(Z,1)),:) ~= Z, 2)-sum(Z==5, 2)-sum(germline.seq==5, 2);
cons_to_germ = sum(cons ~= germline.seq)-sum(cons == 5)-sum(germline.seq == 5)

thresh = 9;
if length(samples) >=2
	cons2 = Z(find(stamps==M+1,1),:);
	dist_to_cons2 = sum(cons2(ones(1,size(Z,1)),:) ~= Z, 2)-sum(Z==5, 2)-sum(cons2==5, 2);
	dist_to_cons2(1)
	cons2_to_germ = dist_to_germ(find(stamps==M+1,1))
	ix = (dist_to_cons+dist_to_germ <= cons_to_germ+2*thresh)  | ...
	(dist_to_cons2+dist_to_germ <= cons_to_germ+2*thresh) | ...
	(dist_to_cons+dist_to_cons2 <= dist_to_cons2(1)+2*thresh);
else
	ix = (dist_to_cons+dist_to_germ <= cons_to_germ+2*thresh);
end

Z = Z(ix,:);
stamps = stamps(ix);

end

function unittest()

samples = [72 73];
cd ~/JVL/src/VDJ
addpath('../DPtrees/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));

[Z stamps germline] = parse_dna(samples);

%%  Run DPTree
    params = [];
    params.me = 0.5;
    params.alpha = 0.5;
    params.epsilon = 0.005;
    params.sigma = 0.025;
    params.germline = germline.seq;    
    params.greedy = 0;
    [T seqs t trace best dbg]= DPtrees_inference_binary(Z, 2000, params);

params = [];
params.me = best.me;
params.alpha = best.alpha;
params.epsilon = best.epsilon;
params.sigma = best.sigma;
params.germline = germline.seq;
params.greedy = -1;
[T_ seqs_ t_ trace_ best_ dbg_] = DPtrees_inference_binary(Z, 50, params, best.T, best.sequences, best.t);

save(sprintf('S%snormalized', sprintf('%d_', samples)), 'best_', 'Z', 'trace', 'trace_', 'best', 'germline', 'stamps', 'samples');

%load S58_normalized; 
h = view_tree(best_.T, best_.sequences, 1, [best_.t; stamps]'); saveas(h.hgAxes, sprintf('S%stree.jpg', sprintf('%d_', samples)));

params = [];
params.me = best_.me;
params.alpha = best_.alpha;
params.epsilon = best_.epsilon;
params.sigma = best_.sigma;
params.germline = germline.seq;


end



function [Z stamps]= parse_dna_file(filename, remove_head, thresh)
% first line is the header
if nargin <2 , remove_head = false; end
if nargin < 3, thresh = 1000; end

X = importdata(filename);

counts = X.data;
reads = X.textdata(2:end, 1);
label = X.textdata{1};


map('ACGTN') = 1:5;
Y = map(cell2mat(reads));
ix = counts(:,1) < thresh;  % ignore reads with distance >=thresh
counts = counts(ix,2:end); % remove first column (distance)
Y = Y(ix,:);

if remove_head % the max count is replaced by just one count
    s = sum(counts, 2)-max(counts, [], 2)+1;
    counts = s;
end

spread = [];
stamps = [];
for i=1:size(Y,1)
    for j=1:size(counts,2)
        spread = [spread i*ones(1,counts(i,j))];    
        stamps = [stamps j*ones(1,counts(i,j))];
    end
end

Z = Y(spread, :);


% second line is the consensus


% other lines are individual reads, already filtered to have the same
% length


end
%%
%saveas(h(2), sprintf('tree%s.jpg', get(get(ch(1), 'SelectedObject'), 'tag')));
