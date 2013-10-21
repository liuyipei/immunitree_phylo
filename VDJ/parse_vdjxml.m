function [rep V D J] = parse_vdjxml(filename, V, D, J, bin_dir)
%function [rep V D J] = parse_vdjxml(filename, V, D, J, bin_dir)
% input:
% filename - path to xml file
% bin_dir - directory to put output files
fid = fopen(filename);
fields = {'descr', 'seq', 'v', 'd', 'j', 'ighc'};

stOut = 1; stIn = 2;
state = stOut; 
ct = 0;
X = [];
if nargin<5
    bin_dir = [];
end

if nargin<4
	V = {}; D = {}; J = {};
	rep = zeros(0, 0, 0);
else
	rep = zeros(length(V), length(D), length(J));
end


fod = cell(length(V), length(J));
%if ~isempty(bin_dir)
%    for v = 1:size(fod,1)
%        for j = 1:size(fod,2)
%            fod{v,j} = fopen(sprintf('%s/bin_%d_%d', bin_dir, v, j), 'w');
%        end
%    end
%end

%try
while 1
    str = fgetl(fid); 
    if ~ischar(str), break, end
    if state == stOut
        if strcmp(str, '<ImmuneChain>')
            %
            % About to start processing a record.
            % 
            state = stIn;
            fld = 1;
            ct = ct + 1;
        end
        continue;
    end
    if state == stIn
        if strcmp(str, '</ImmuneChain>')
            state = stOut;
            if mod(ct, 10000) == 0, fprintf('%d\n', ct); end
            assert(fld == length(fields)+1);
            %
            % Just finished processing a record.  X has it.
            % 
            vix = find(strcmp(X.v, V));
            if isempty(vix) %&& ~isempty(X.v)
                V{end+1} = X.v;
                vix = length(V);
                rep(vix, 1, 1) = 0;
            end
            dix = find(strcmp(X.d, D));
            if isempty(dix) %&& ~isempty(X.d)
                D{end+1} = X.d;
                dix = length(D);
        		rep(1, dix, 1) = 0;
            end
            jix = find(strcmp(X.j, J));
            if isempty(jix) %&& ~isempty(X.j)
                J{end+1} = X.j;
                jix = length(J);
        		rep(1, 1, jix) = 0;
            end
            rep(vix, dix, jix) = rep(vix, dix, jix)+1;

            if ~isempty(bin_dir)
                if isempty(fod{vix, jix})
                    fod{vix,jix} = fopen(sprintf('%s/bin_%d_%d.fasta', ...
                                                 bin_dir, vix, jix), 'w');
                end                
                fprintf(fod{vix, jix}, '>%s,%s,%s,%s,%s,%s\n%s\n', ...
                        X.descr, X.d, X.ighc,  X.barcode, X.experiment, X.clone, X.seq);
            end
%            write_to_file(fod{vix, jix}, X);            
            continue;
        end            
        ix1 = find(str == '<');
        ix2 = find(str == '>');
        if (length(ix1) == 2) && (length(ix1) == 2)
            field = str(ix1(1)+1:ix2(1)-1);
            if fld <= length(fields)
                field_ = fields{fld};
                assert(strcmp(field, field_));
                val = str(ix2(1)+1:ix1(2)-1);
                X = setfield(X, field, val);
                fld = fld + 1;
            else
                if strcmp(field, 'tag')
                    val = str(ix2(1)+1:ix1(2)-1);
                    ix = find(val=='|', 1);
                    if ~isempty(ix) && strcmp(val(1:ix(1)-1), 'experiment')
                        X.experiment = val(ix(1)+1:end);
                    end
                    if ~isempty(ix) && strcmp(val(1:ix(1)-1), 'clone')
                        X.clone = val(ix(1)+1:end);
                    end                                   
                    if ~isempty(ix) && strcmp(val(1:ix(1)-1), 'barcode')
                        X.barcode = val(ix(1)+1:end);
                    end
                end
            end
        end
        continue;
    end        
    
end

%catch exception   
%    fprintf('Exception\n');
%    fclose(fid);
%end
   
fclose(fid);
if ~isempty(bin_dir)
    for v = 1:size(fod,1)
        for j = 1:size(fod,2)
            if ~isempty(fod{v,j})
                fclose(fod{v,j});
            end
        end
    end
end
    
end


%function write_to_file(fid, X)
%end

