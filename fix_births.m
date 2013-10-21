function tree = fix_births(tree, F)
    b = [tree(:,2); F];
    [~,ord] = sort(b);
    for i=2:(length(ord)-1)
        if b(ord(i))-b(ord(i-1)) < 1e-9
            if b(ord(i+1))-b(ord(i)) > 1e-9
               b(ord(i)) = b(ord(i))+1e-9; 
            else
               b(ord(i)) = (b(ord(i-1))+b(ord(i+1)))/2;
            end
        end
    end
    tree(:,2) = b(1:end-1);
end


function test()
%%
%fix_births(tree_,10)-tree_
junk = [0 1 0; 0 2 0; 0 1+1e-10 0; 0 2 0; 0 2+1e-10 0; 0 3 0; 0 1.5 0; 0 2+5e-10 0];
%junk_ = diff(sort(junk(:,2)))
junk_ = fix_births(junk,3);
[~,ord] = sortrows(junk);
for i=1:length(junk), fprintf('%.11f  %.11f\n', junk(ord(i),2), junk_(ord(i),2)); end
fprintf('\n');

end