function [out1,out2] = get_wavevectors(N_wv,a,options)
    
    if options.isTrimRightBoundary
        wv = get_IBZ_wavevectors(N_wv + [1 0],a,'none',1); % Add one pixel of resolution to x direction since it will get trimmed off later.
%         wv = get_IBZ_wavevectors(N_wv,a,'none',1); % *DON'T* add one pixel of resolution to x direction even though one will get trimmed off later.
        wv(wv(:,1) == pi/a,:) = [];
    elseif ~options.isTrimRightBoundary
        wv = get_IBZ_wavevectors(N_wv,a,'none',1);
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

