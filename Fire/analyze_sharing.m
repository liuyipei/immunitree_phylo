%% Get the (igBlast) V-D-J of every reads from the header

vdj_file_name = [fasta_file_name '.igblast'];
fid = fopen(vdj_file_name);
junk = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
%%
[data,~,~,ix] = correct_PCR_artifacts(data);
junk = junk{1};
junk = junk(ix);
[~,I,J] = unique(junk);

%% overlap between reads across different replicates
HIT = zeros(240);
map = zeros(24,1); map([11 12 21 22 23 24]) = 1:6;
for i=1:length(I)
    if mod(i, 10000) == 0 , fprintf('%d, ', i); end
    ix = (J==i);
    replicas = 6*(double(data.subject_num(ix))-1) + map([double(data.FR(ix)*10) + data.rep(ix)-'a'+1]);
    replicas = unique(replicas); 
    for k=1:length(replicas)
        for l=k+1:length(replicas)
            HIT(replicas(l),replicas(k)) = HIT(replicas(l),replicas(k)) + 1;
        end
    end    
end
fprintf('\n');
imagesc(HIT, [0 20]);
plot_class_lines(1:6:241, false, 'k', true);
plot_class_lines(1:6:241, true, 'k', true);
plot_class_lines(1:24:241, false, 'm', false);
plot_class_lines(1:24:241, true, 'm', false);
%% overlap between reads across different samples
HIT2 = zeros(40);
for i=1:length(I)
    if mod(i, 10000) == 0 , fprintf('%d, ', i); end
    ix = (J==i);
    replicas = unique(data.subject_num(ix)); 
    for k=1:length(replicas)
        for l=k+1:length(replicas)
            HIT2(replicas(l),replicas(k)) = HIT2(replicas(l),replicas(k)) + 1;
        end
    end    
end
fprintf('\n');
%%
figure;
imagesc(HIT2, [0 200]);
plot_class_lines(1:4:241, false, 'm', true);
plot_class_lines(1:4:241, true, 'm', true);
colorbar('east');