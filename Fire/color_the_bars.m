function color_the_bars(h)
    if length(h) == 4
        set(h(1),'FaceColor',[1 0 0]);
        set(h(2),'FaceColor',[1 1 0]);
        set(h(3),'FaceColor',[0 1 1]);
        set(h(4),'FaceColor',[0 0 1]);
    end
    if length(h) == 6
        jet6 = flipud(jet(6));
        for i=1:length(h), 
            set(h(i), 'FaceColor', jet6(i,:));
        end    
    end
end