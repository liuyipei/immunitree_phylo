function a = add_germline_as_root(a, germline)
    a.tree(:,1) = a.tree(:,1)+1;
    a.tree = [ 0 0 a.tree(1,3)  ;a.tree];
    a.sequences = [germline ; a.sequences];
    a.t = a.t+1;
    a.tree2(:,1) = a.tree2(:,1)+1;
    a.tree2 = [ 0 0 a.tree2(1,3)  ;a.tree2];

end