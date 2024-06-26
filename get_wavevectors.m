function [out1,out2] = get_wavevectors(N_wv,a,options)
    
    if options.isTrimRightBoundary & strcmp(options.format,'grid')
        error('Cannot simultaneously trim right boundary AND output in grid format.')
    end
    
    wv = get_IBZ_wavevectors(N_wv,a,'none',1); % Add one pixel of resolution to x direction since it will get trimmed off later.
    if options.isTrimRightBoundary
        wv(wv(:,1) == pi/a,:) = [];
        wv(wv(:,2) == pi/a & wv(:,1) > 0,:) = [];
        wv(wv(:,2) == 0 & wv(:,1) > 0,:) = [];
    elseif ~options.isTrimRightBoundary
        % Do nothing
    end
    
    if strcmp(options.format,'grid')
        WV = cat(3,reshape(wv(:,1),[N_wv(2) N_wv(1)]),reshape(wv(:,2),[N_wv(2) N_wv(1)]));
        out1 = WV(:,:,1);
        out2 = WV(:,:,2);
    elseif strcmp(options.format,'list')
        % do nothing, it's already a list
        out1 = wv;
    end
end

