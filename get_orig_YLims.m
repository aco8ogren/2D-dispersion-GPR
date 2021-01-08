function orig_YLims = get_orig_YLims(ax)
    [m,n] = size(ax);
    for i = 1:m
        for j = 1:n
            orig_YLims(:,i,j) = ax(i,j).YLim;
        end
    end
end