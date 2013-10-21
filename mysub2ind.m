%function ix = mysub2ind(B, x, y, z)
% same as sub2ind([B B L], x, y, z)
% put z=1 for L=1.
function ix = mysub2ind(B, x, y, z)
    % have to convert to double! Otherwise limited by int type
    x = double(x); y = double(y); z = double(z);
    ix = x + (y-1)*B + (z-1)*B^2;
end
% 
% function ix = mysub2ind(B,x,y)
%     ix = x + (y-1)*B;
% end


function test()
%%
B = 64;
L = 100;
N = 30;
x = int16(ceil(B*rand(1,N)));
y = int16(ceil(B*rand(1,N)));
z = int16(ceil(L*rand(1,N)));
ix = mysub2ind(B, x, y, z);
assert(isequal(sub2ind([B B L], double(x), double(y), double(z)), ix));
fprintf('Test Passed\n');

end