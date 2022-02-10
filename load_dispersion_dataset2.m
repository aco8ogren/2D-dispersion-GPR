function data = load_dispersion_dataset2(data_path)
    data = load(data_path);
    
    wv = data.WAVEVECTOR_DATA(:,:,1);
    [~,idxs] = sort(wv(:,2));
    data.WAVEVECTOR_DATA = data.WAVEVECTOR_DATA(idxs,:,:);
    data.EIGENVALUE_DATA = data.EIGENVALUE_DATA(idxs,:,:);
    if isfield(data,'EIGENVECTOR_DATA') % Not yet sure if this is right or if i should do ~isempty(data.EIGVECDATA)
        data.EIGENVECTOR_DATA = data.EIGENVECTOR_DATA(idxs,:,:,:);
    end
end