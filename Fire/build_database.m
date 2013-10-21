%% Analyze a new read file
fasta_file_name = '/afs/cs/u/joni/scratch/data/lymph/t5.fa';
Z = fastaread(fasta_file_name);
Z_ = struct2cell(Z);
data.reads = Z_(2,:)';

%%  process data
N = length(Z);
str = cell2mat(cellfun(@(x) x([1:5 7 8 15 16 18:25]), Z_(1,:)', 'uniformoutput', false));
%str = cell2mat(cellfun(@(x) ['5_' x([1:3 5 6 13 26])], Z_(1,:)', 'uniformoutput', false));
[ids I J] = unique(str(:,1:5), 'rows');
h = hist(J, 1:max(J)); % = diff([0 I'])
data.RF = str2num(str(:,6:7)); 
%%
%if (length(unique(RF)) == 1), clear RF; else assert(false); end
data.FR = str(:,8)-'0';
data.rep = str(:,9);
data.subjects = ids;
data.subject_num = J;
data.read_id = str(:,10:17);

%% this works for tibshirani and utsw
% Use textread or any other function to read the header mapping files
% tibheadermap.txt
% utswheadermap.txt
metadata_file_name = [fasta_file_name '.map'];
fid = fopen(metadata_file_name);
junk = textscan(fid, '%s %s %s %s %s %s', 'delimiter', '\t' );
if 1, junk = junk([1 6 5 2 3 4]); end  % the .map file for t5 is different from the tibshirani file.
fclose(fid);
%% get meta data and statistics on subjects
subjects = struct('code',[], 'id', [], 'sample', [], 'source', [], 'desc',[]);
subjects = subjects(ones(size(ids,1),1));
for i=1:size(ids,1)
    i
    ix = ~cellfun(@isempty, strfind(junk{1}, ids(i,:)));       
    iy = (J == i);   

    subjects(i).code = ids(i,:);    
    
    res = unique(junk{2}(ix,:));
    assert(size(res,1) == 1);
    subjects(i).sample = res{1};

    res = unique(junk{4}(ix,:));
    assert(size(res,1) == 1);
    subjects(i).id = res{1};
    
    res = unique(junk{5}(ix,:));
    assert(size(res,1) == 1);
    subjects(i).source = res{1};

    res = unique(junk{6}(ix,:));
    assert(size(res,1) == 1);
    subjects(i).desc = res{1};

    subjects(i).reads = [sum(ix) h(i)];
       
end
data.subjects = subjects;
%save([fasta_file_name '.mat'], 'data');
%save([fasta_file_name '.mat'], 'data', 'subjects');

%% Get the (igBlast) V-D-J of every reads from the header
global rep
rep = load_repertoire('igblast');

vdj_file_name = [fasta_file_name '.igblast'];
fid = fopen(vdj_file_name);
junk = textscan(fid, '%s %s %s %s %s', 'delimiter', '_' );
fclose(fid);
igblast.V = junk{1};
igblast.D = junk{2};
igblast.J = junk{3};

data.igblast = zeros(length(data.reads),3);
[TF, data.igblast(:,1)] = ismember(igblast.V, {rep.V.Header});
[TF, data.igblast(:,2)] = ismember(igblast.D, {rep.D.Header});
[TF, data.igblast(:,3)] = ismember(igblast.J, {rep.J.Header});


%% Get the (ihmmune) V-D-J of every reads from the header
global rep
rep = load_repertoire('ihmmune');
vdj_file_name = [fasta_file_name '.ihmmune'];
fid = fopen(vdj_file_name);
junk = textscan(fid, '%s %s %s %s', 'delimiter', '\t' );
fclose(fid);
ihmmune.id = junk{1};
ihmmune.V = junk{2};
ihmmune.D = junk{3};
ihmmune.J = junk{4};
%% map records to reads
Z = fastaread([fasta_file_name '.headers']);
[TF, ids] = ismember(ihmmune.id, {Z.Header}); % 10 secs

%% collect all IGH* genes
for c = 'VDJ';
    tic
    loc = zeros(100, 2);
    str = cell(100,1);
    k = 0;
    for i=1:length(ihmmune.(c))
        t = 0;    
        x = ihmmune.(c){i};
        ix = strfind(x, sprintf('IGH%c', c));
        if k > length(str)-10
            str{2*length(str)} = [];
            loc(2*length(str),2) = 0;
            fprintf('%d records\n', k);    
            toc
        end
        for j=ix    
            if j>1 && x(j-1)~= ' ', continue; end
            y = strtok(x(j:end), '] \t');
            t = t+1;
            k = k+1;
            loc(k,:) = [ids(i) t];
            str{k} = y;
        end
    end
    % Using the given reprtoire
    rep = load_repertoire('ihmmune');
    [TF, v] = ismember(str(1:k), {rep.(c).Header}); 
    ihmmune_inferred.(c) = full(sparse(loc(1:k,1), loc(1:k,2), v(1:k), length(Z),max(loc(:,2))));
    sum(TF == 0)
end
ihmmune = ihmmune_inferred;
clear ihmmune_inferred
%%  ihmmune: Save data and clean V headers from duplicates
save([fasta_file_name '.ihmmune.mat'], 'rep', 'ihmmune');
data.ihmmune = [ihmmune.V(:,1) ihmmune.D ihmmune.J ihmmune.V(:,2:end)];

% clean duplicate V entries for same read
ix = []; 
y = data.ihmmune(:,[1 4:7]); 
for i=1:size(y, 1) 
    if mod(i,10000) == 0, fprintf('%d, ', i); end; 
    x = y(i, y(i,:) > 0);  
    x_ = unique(x);
    if length(x_) ~= length(x)        
        x_ = [x_ zeros(1,size(y,2)-length(x_))];
        y(i,:) = x_;
        ix = [ix i]; 
    end; 
end
fprintf('\n corrected %d headers with duplicates.\n', length(ix));
data.ihmmune(:, [1 4:7]) = y;
save([fasta_file_name '.mat'], 'data');



%% Creating the repertoire from scratch
% no need to run this segment as I got the correct catalogue from Katherine
rep.code = 'ihmmune_inferred';
N = length(Z);
[rep.(c), ~, v] = unique(str(1:k));
ihmmune_inferred.(c) = full(sparse(loc(1:k,1), loc(1:k,2), v(1:k), N,max(loc(:,2))));


%% compute the V-D-J of every read
global rep
rep = load_repertoire('igblast');
map('ACGTN') = 1:5;
N = length(data.reads);
vdj = zeros(N,3);
progress = ceil(N/100);
for i=1:N
    if mod(i,progress) == 0, fprintf('*'); end
    x = map(data.reads{i});
    try
        [vdj(i,1),vdj(i,3),~, ~, vdj(i,2)] = analyze_VDJ(x);
    end
end
data.vdj = vdj;
save([fasta_file_name '.mat'], 'data');




%  %  Done making .mat data file!  %  %


%%  Find redundant V ids
N = size(ihmmune.V,1);
M = size(rep.V,1);
I = false(M,M);
ix = false(N,M);
for i=1:M
    ix(:,i) = sum(ihmmune.V == i, 2) > 0;
end
for i=1:M
    if mod(i,10)==0, i, end
    for j=i+1:M
        I(i,j) = all(ix(:,i) <= ix(:,j)); % every time j is listed so is i
        I(j,i) = all(ix(:,j) <= ix(:,i)); % every time i is listed so is j
    end
end

redundant = find(sum(I,2)' > 0 & sum(ix) > 0)



