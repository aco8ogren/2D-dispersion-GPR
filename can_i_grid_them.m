clear; close all;

wv_i_x = linspace(-1,1,5)';
wv_i_y = linspace(0,1,3)';

wv_j_x = linspace(-1,1,11)';
wv_j_y = linspace(0,1,6)';

[WV_i_x,WV_i_y] = meshgrid(wv_i_x,wv_i_y);

[WV_j_x,WV_j_y] = meshgrid(wv_j_x,wv_j_y);

[ndWV_i_x,ndWV_i_y,ndWV_j_x,ndWV_j_y] = ndgrid(wv_i_x,wv_i_y,wv_j_x,wv_j_y);

disp('wv_i')
size(WV_i_x)
prod(size(WV_i_x))

disp('wv_j')
size(WV_j_x)
prod(size(WV_j_x))

disp('nd')
size(ndWV_i_x)
prod(size(ndWV_i_x))

counter = 0;
for i = 1:size(WV_i_x,1)
    for j = 1:size(WV_i_x,2)
        counter = counter + 1;
        wv_i_c(counter,:) = [WV_i_x(i,j) WV_i_y(i,j)];
    end
end

counter = 0;
for i = 1:size(WV_j_x,1)
    for j = 1:size(WV_j_x,2)
        counter = counter + 1;
        wv_j_c(counter,:) = [WV_j_x(i,j) WV_j_y(i,j)];
    end
end

WV_c = combvec(wv_i_c',wv_j_c');

disp('using combvec')
size(WV_c)

% I think I can grid them

% wv_q = combvec(wv_i',wv_j')