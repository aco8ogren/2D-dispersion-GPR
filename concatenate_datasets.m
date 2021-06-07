clear; close all;

fns{1} = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\ground_truth output 20-May-2021 17-11-07\DATA N_struct128 N_k201 RNG_offset0 20-May-2021 17-11-07.mat";
fns{2} = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\2D-dispersion\OUTPUT\ground_truth2 output 21-May-2021 10-32-17\DATA N_struct100 N_k201 RNG_offset0 21-May-2021 10-32-17.mat";

N_dsets = length(fns);

key_cell = {'modulus','density','poisson'};
CONSTITUTIVE_DATAS = containers.Map(key_cell,{cell(0,0),cell(0,0),cell(0,0)});
for key_idx = 1:length(key_cell)
    key = char(key_cell{key_idx});
    CONSTITUTIVE_DATAS(key) = cell(N_dsets,1);
end
EIGENVALUE_DATAS = cell(N_dsets,1);
WAVEVECTOR_DATAS = cell(N_dsets,1);

for fn_idx = 1:N_dsets
    fn = fns{fn_idx};
    load(fn);
    EIGENVALUE_DATAS{fn_idx} = EIGENVALUE_DATA;
    WAVEVECTOR_DATAS{fn_idx} = WAVEVECTOR_DATA;
    for key_idx = 1:length(key_cell)
        key = key_cell{key_idx};
        temp{fn_idx,key_idx} = CONSTITUTIVE_DATA(key);
    end
end

EIGENVALUE_DATA = cat(3,EIGENVALUE_DATAS{:});
WAVEVECTOR_DATA = cat(3,WAVEVECTOR_DATAS{:});

% CONSTITUTIVE_DATA = containers.Map(key_cell,{zeros(0),zeros(0),zeros(0)});
% for key_idx = 1:length(key_cell)
%     key = key_cell{key_idx};
% %     CONSTITUTIVE_DATAS(key) = temp{:,key_idx};
% %     temp2 = CONSTITUTIVE_DATAS(key);
%     CONSTITUTIVE_DATA(key) = cat(3,temp{:,key_idx});
% end

CONSTITUTIVE_DATA = 'couldn''t concatenate because designs are different resolutions';

save('concatenated_dataset','EIGENVALUE_DATA','WAVEVECTOR_DATA','CONSTITUTIVE_DATA')

text_to_write_str = cat(1,"This dataset is a concatenation of ",fns{:});
text_to_write = [];
for i = 1:size(text_to_write_str,1)
    text_to_write = [text_to_write regexprep(char(text_to_write_str(i,:)),'\','\\\') ' ' newline];
end
text_to_write = strtrim(text_to_write);

out_fn = 'concatenation_details.txt';
fileID = fopen(out_fn,'w');
fprintf(fileID,text_to_write);
fclose(fileID);
fclose('all');



