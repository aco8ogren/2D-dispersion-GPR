function make_y_axis_uniform(ax,max_dir)
    [m,n] = size(ax);
    for i = 1:m
        for j = 1:n
            YLims(:,i,j) = ax(i,j).YLim;
        end
    end
    
    YMins = min(YLims(1,:,:),[],max_dir);
    YMaxs = max(YLims(2,:,:),[],max_dir);
    
    if max_dir == 1;
        for i = 1:m
            for j = 1:n
                set(ax(i,j),'YLim',[YMins(i) YMaxs(i)])
            end
        end
    elseif max_dir == 2;
        for i = 1:m
            for j = 1:n
                set(ax(i,j),'YLim',[YMins(j) YMaxs(j)])
            end
        end
    elseif strcmp(max_dir,'all')
        for i = 1:m
            for j = 1:n
                set(ax(i,j),'YLim',[YMins YMaxs])
            end
        end
    end
    
end