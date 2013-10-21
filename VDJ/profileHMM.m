addpath('../DPtrees/');
addpath(genpath('/afs/cs/u/joni/scratch/software/lightspeed'));

%%
% choose a V D J
rep.V = fastaread('V.fa');
rep.D = fastaread('D.fa');
rep.J = fastaread('J.fa');
% choose V, D, J.
v = ceil(length(rep.V)*rand);
d = ceil(length(rep.D)*rand);
j = ceil(length(rep.J)*rand);
V_seq = rep.V(v).Sequence; 
D_seq = rep.D(d).Sequence; 
J_seq = rep.J(j).Sequence;

fprintf('Generating sequence...\n');
[seq V N1 D N2 J ] = profileHMM_gen(V_seq, D_seq, J_seq);
fprintf('Annotating sequence...\n');
[V_ N1_ D_ N2_ J_] = profileHMM_annotate(seq, V_seq, D_seq, J_seq);

fprintf('V\n%s\n%s\n%s\n', V_seq, V, V_);
fprintf('N1\n%s\n%s\n', N1, N1_);
fprintf('D\n%s\n%s\n%s\n',D_seq, D, D_);
fprintf('N2\n%s\n%s\n', N2, N2_);

str_ = ['%' num2str(length(J_seq)) 's'];
str = sprintf('J\n%s\n%s\n%s\n', str_, str_, str_);
fprintf(str, J_seq, J, J_);

%%
me = 0.5; alpha = 0.5; epsilon = 0.01; sigma = 0.02; 
[reads, sequences, T, nodes] = generate_data(100, seq, me, alpha, epsilon, sigma);
%%
DPtrees_inference_binary(reads, 0, me, alpha, T, sequences, nodes);
%%
[T_ sequences_ t_ trace] = DPtrees_inference_binary(reads, 100, me, alpha);
% start creating a mutation tree
%%  P_del stats
load('P_del_stats');
stats = {V_stats, D1_stats, D2_stats, J_stats};
names = {'V', 'D1', 'D2', 'J'};
for i=1:4
    P_del(i).title = names{i};
    P_del(i).pdf_ex = hist(stats{i},0:20)/length(stats{i});
    P_del(i).cdf_ex = [fliplr(cumsum(fliplr(P_del(i).pdf_ex))) 0]; % Prob T >= i-1
    P_del(i).inc_ex = P_del(i).cdf_ex(2:end)./P_del(i).cdf_ex(1:end-1);    
    % plot
    subplot(2,2,i);
    bar(0:20, P_del(i).pdf_ex);
    title(P_del(i).title);
end

% %% N_add stats
% figure; 
% subplot(1,2,1);
% plot(0:40, N_add(1).pdf);
% title('N1');
% ylim([0 0.15]);
% subplot(1,2,2);
% plot(0:40, N_add(2).pdf);
% ylim([0 0.15]);
% title('N2');
%%
tic
hmms = cell(nSeq, 1);
for s=1:nSeq
    hmm{s} = hmmprofstruct(length(seqs_nt{s}), 'Alphabet', 'NT');
    hmm{s}.MatchEmission = exp(p(seqs_nt{s}, 1:4));        
end
for i=1:N, i, x = int2nt(X{i}); for s=1:nSeq, hmmprofalign(hmm{s}, x); end; end; toc

%%
tic
for i=1:N, x = int2nt(X{i}); for s=1:nSeq, nwalign(seqs{s}, x); end; end; toc

