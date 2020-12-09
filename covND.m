function out = covND(A)   
    s = size(A);
    dim_of_data = length(s) - 1;       
    out = zeros(s(1)*ones(dim_of_data,1));   
    
    means = mean(A,dim_of_data + 1);
    for s_idx = 1:dim_of_data
        for idx = 1:s(s_idx)
            
end