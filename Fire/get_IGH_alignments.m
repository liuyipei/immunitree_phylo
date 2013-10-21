function A = get_IGH_alignments(rep, c)

Z = fastaread(sprintf('IGH%c.fa', c));
assert(length(Z) == length(rep.(c)));
l = cellfun(@length, {Z.Sequence});
A = false(length(Z), max(l));

fprintf('\nmismatches between alignment file the %c repertoire:\n', c);
for i=1:length(Z), 
    x = char(Z(i).Sequence(Z(i).Sequence~='.'));  % IMGT
    if x(1) >= 'a'
        x = x+'A'-'a';  % convert to upper case (if needed)
    end
    y = rep.(c)(i).Sequence; % our repertoire
    if ~isequal(x, y)
        l = min(length(x),length(y));
        delta = sum(x(1:l)~=y(1:l));
        delta_len = length(x)-length(y);
        fprintf('%d %s %s\tdiff = %d\tlen_diff = %d\n', i,Z(i).Header(1:min(end,22)), rep.(c)(i).Header, delta, delta_len); 
        assert(abs(delta_len)<3);
        assert(delta<3);
        if delta_len > 0 % y (repertoire) is longer
            assert( all(Z(i).Sequence(end-delta_len+1 : end)~= '.'));
            Z(i).Sequence = Z(i).Sequence(1:end-delta_len);
        elseif delta_len < 0 % x (IMGT) is longer
            Z(i).Sequence = [Z(i).Sequence rep.(c)(i).Sequence(end+delta_len+1 : end)];
        end        
    end
%     rep.V(i).al = false(1,max_l);
%     rep.V(i).al(Z(i).Sequence ~= '.') = true;
    A(i, Z(i).Sequence ~= '.') = true;
end
fprintf('Good enough!\n');

end