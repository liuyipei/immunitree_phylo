function [nseq dmap] = apply_edit(seq,ed)
%function nseq = apply_edit(seq,ed)
%  0 - match
%  1 - insert
% -1 - delete
% Also returns a del map.  an length(seq)+1 array.  for each '-1' it finds the
% nearest '0' on the right and puts a '1' there.

    nseq = zeros(1,length(ed));
    dmap = zeros(1,length(ed)+1);;
    ct = 1; lct = 1;
    del_count = 0;
    for i=1:length(ed)

        if ed(i) == 0
            nseq(lct) = seq(ct);
            dmap(lct) = del_count;
            ct = ct+1;
            lct = lct+1;
            del_count = 0;
        else
            if ed(i) == -1
                del_count = del_count + 1;
                ct = ct+1;
            else %  if ed(i) > 0
                nseq(lct) = 5;%ed(i);
                lct = lct+1;
            end            
        end
    end
    nseq = nseq(1:lct-1);
    dmap = dmap(1:lct);
    dmap(lct) = del_count;
end

