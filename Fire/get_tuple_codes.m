% index which combinations of sources are "pure", "doubles" and "triples"
function [codes singles doubles triplets] = get_tuple_codes()
    codes = (dec2bin(0:15)=='1');
    [~, singles] = ismember(eye(4), codes,'rows');
    doubles = flipud(find(sum(codes,2) == 2));
    [~, triplets] = ismember(1-eye(4), codes, 'rows');
end