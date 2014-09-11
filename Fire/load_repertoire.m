function rep = load_repertoire(code, replace_underscores, VDJ_letters)
    
    if nargin<2, replace_underscores = false; end
    if nargin<3, VDJ_letters = 'VDJ'; end

    % [~,whoami] = system('whoami')
    liuyipei = 1; %isequal(whoami(1:min(8,end)), 'liuyipei')
    %rep_dir = '/afs/cs/u/joni/scratch/data/VDJrepertoire/';
    if liuyipei
        rep_dir = '.'; %
    end
    for c = VDJ_letters
        sprintf([rep_dir '/' code '/%c.fa'], c)
        rep.(c) = fastaread(sprintf([rep_dir '/' code '/%c.fa'], c));
        % removing spaces from headers
        for i=1:length(rep.(c))
            rep.(c)(i).Header = strtok(rep.(c)(i).Header, ' ');

            % replace '_' with '-' so that '_' can be used as delimeter
            % also replace '-' with '*'.
            if replace_underscores  
                fprintf('Replacing underscored...\n');
%                rep.(c)(i).Header(rep.(c)(i).Header == '_') = '*';
                str = rep.(c)(i).Header;
                if isempty(find(str == '*',1))
                    assert(sum(str == '-') <= 1);
                    str(str=='-') = '*';
                end
                str(str=='_') = '-';
                rep.(c)(i).Header = str;
            end

            % convert to capitals and remove all letters that are not ACGTN
            tmp = rep.(c)(i).Sequence;  
            tmp(tmp>='a') = tmp(tmp>='a') -'a'+'A'; 
            ix = ~ismember(tmp, 'ACGTN');
            tmp(ix) = 'N';
            rep.(c)(i).Sequence = tmp;
        end
    end    
    rep.code = code;
    
    

end
