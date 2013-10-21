%%  Add to path
cd /vision/u/liuyipei/joni/phylo/Fire
phylo_path = '/vision/u/liuyipei/joni/phylo/';
addpath(phylo_path, '-end');
addpath([phylo_path 'Fire'], '-end'); % fixes issues with affinegapmex/nwalign
addpath([phylo_path 'VDJ'], '-end');
addpath([phylo_path 'util'], '-end');

%%
%"pt45 HC"   -Patient 45 Heavy Chain Variable(SRR275668);
%"pt45 kappa"-Patient 45 Light Chain Variable(SRR275679);
%"pt74 HC set1"-Patient 74 Heavy Chain Variable(SRR275711);
%"pt74 HC set2"-Patient 74 Heavy Chain Variable(SRR277211);
hiv_dir = '/vision/u/liuyipei/vdj_data/wu2011science/';
files = {'SRR275668.fasta', 'SRR275679.fasta', 'SRR275711.fasta', 'SRR277211.fasta'};
i = 1
file = files{1}
chain_type = 1
if i == 3
  chain_type = 0
end
[a h chain] = pipeline_for_hiv([hiv_dir file], [], chain_type, true)
