function set_YLims(ax,YLims)
    [m,n] = size(ax);
    for i = 1:m
        for j = 1:n
            set(ax(i,j),'YLim',YLims(:,i,j));
        end
    end
end