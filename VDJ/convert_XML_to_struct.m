function X = convert_XML_to_struct(st)
    
%st = parseXML('/afs/cs/u/joni/scratch/data/Uri/delme.xml');
%st.Children = st.Children(2:2:end);
fields = {st.Children(2).Children(2:2:end).Name};
assert( strcmp(fields{10}, 'tag'));
assert(~strcmp(fields{9},  'tag'));
X = [];
for j=1:9
    X = setfield(X, fields{j}, []);
end
    X = setfield(X, 'tag', {});

ct = 0;
for i=2:2:length(st.Children)
    ct = ct + 1;
    for j=2:2:length(st.Children(i).Children)
        field = st.Children(i).Children(j).Name;
        if ~isempty(st.Children(i).Children(j).Children)
            val = st.Children(i).Children(j).Children(1).Data;
            if ~strcmp(field, 'tag')
                X = setfield(X,{ct},field, val);
            else
                    X(ct).tag = [X(ct).tag {val}];
            end
        end
    end
end

end