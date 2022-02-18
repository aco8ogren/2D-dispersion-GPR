classdef EmpiricalCovarianceFunctionConstructor
    properties
        dataset_path
        default_query_format char
        precomputed_combvec
        covariance_matrix
        grid_vectors
        wavevector_array_size
        gridded_interpolant
        interpolation_method char;
    end
    methods
        function obj = EmpiricalCovarianceFunctionConstructor(dataset_path)
            obj.dataset_path = dataset_path;
        end
        function obj = run(obj)
            data = load(obj.dataset_path);
            dispersion_dataset = data.dispersion_dataset;
            obj.wavevector_array_size = dispersion_dataset.wavevector_array_size;
            obj.grid_vectors = {sort(unique(dispersion_dataset.wavevectors(:,1))),sort(unique(dispersion_dataset.wavevectors(:,2)))};
            for band_index = 1:dispersion_dataset.N_band
                obj.covariance_matrix{band_index} = cov(dispersion_dataset.frequency(:,:,band_index));
            end
            if isempty(obj.interpolation_method)
                obj.interpolation_method = 'cubic';
            end
            for band_index = 1:dispersion_dataset.N_band
                obj.gridded_interpolant{band_index} = griddedInterpolant(obj.grid_vectors,obj.covariance_matrix{band_index},obj.interpolation_method);
            end
        end
        function C = evaluate(obj,wv_i,wv_j,query_format)
            if nargin == 3 % nargin counts obj as an input, so this condition is true if query_format is not specified
                query_format = obj.default_query_format;
            end
            switch query_format
                case 'gridded'
                    % Requires wavevectors to be input sorted by wv_x (for both wv_i and wv_j)
                    wv_i_x = sort(unique(wv_i(:,1)));
                    wv_i_y = sort(unique(wv_i(:,2)));
                    wv_j_x = sort(unique(wv_j(:,1)));
                    wv_j_y = sort(unique(wv_j(:,2)));
                    C4D = C_griddedInterpolant({wv_i_x,wv_i_y,wv_j_x,wv_j_y}); 
                    C4D = permute(C4D,[2 1 4 3]);
                    C = reshape(C4D,size(wv_i_x,1)*size(wv_i_y,1),size(wv_j_x,1)*size(wv_j_y,1));
                case 'scattered'
                    wv = combvec2(wv_i',wv_j'); % gives a 4 x N_combinations array. Faster, but less general than MatLab's combvec.
                    C4D = C_griddedInterpolant(wv');
                    [N_h,~] = size(wv_i);
                    [M_h,~] = size(wv_j);
                    C = reshape(C4D,N_h,M_h);
                case 'scattered - precomputed_combvec'
                    wv = wv_i.values;
                    N_h = wv_i.size(1);
                    M_h = wv_i.size(2);
                    C4D = C_griddedInterpolant(wv');
                    C = reshape(C4D,N_h,M_h);
                otherwise
                    error('query_format not recognized')
            end
        end
    end
end