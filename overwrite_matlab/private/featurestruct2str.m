function str = featurestruct2str(f,feat,width)
%FEATURESTRUCT2STR

% Copyright 2006-2009 The MathWorks, Inc.

fnames = fieldnames(f);
cstr = cell(100,1);

cstr(1) = {[feat blanks(15) f.Location]};
j = 1;
for i = 3:numel(fnames) %start at 3 to save next if
   %if ~any(strcmp(fnames(i),{'Location','Indices'})) 
       if ~iscell(f.(fnames{i}))
           f.(fnames{i}) = {f.(fnames{i})};
       end
       for k = 1:numel(f.(fnames{i}))
           txt = f.(fnames{i}){k};
           if ~isempty(txt) && ischar(txt)
              j = j + 1;
              cstr{j} = ['  /' fnames{i} ' =  ' txt];         
              while numel(cstr{j})>width
                  if strcmp(fnames{i},'translation')
                      cstr{j} = [cstr{j}(1:width-3),'...'];
                  else
                      j = j + 1;
                      sp = find(isspace(cstr{j-1}(1:min(width,end))),1,'last');
                      if isempty(sp) || sp<=10
                          cstr{j} = [blanks(12) cstr{j-1}(width+1:end)];
                          cstr{j-1} = cstr{j-1}(1:width-1);
                      else
                          cstr{j} = [blanks(12) cstr{j-1}(sp+1:end)];
                          cstr{j-1} = cstr{j-1}(1:sp-1);
                      end
                  end
              end
           end
       end
   %end
end
str = char(cstr(~cellfun('isempty',cstr)));

