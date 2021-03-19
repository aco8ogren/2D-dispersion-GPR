function [WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path)
    data = load(data_path,'EIGENVALUE_DATA','WAVEVECTOR_DATA');
    
    EIGENVALUE_DATA = data.EIGENVALUE_DATA;
    WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;
    
    wv = WAVEVECTOR_DATA(:,:,1);
    [~,idxs] = sort(wv(:,2));
    WAVEVECTOR_DATA = WAVEVECTOR_DATA(idxs,:,:);
    EIGENVALUE_DATA = EIGENVALUE_DATA(idxs,:,:);
end