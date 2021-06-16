clear;
N_wv = [52 27];
wv = get_wavevectors(N_wv,1,struct('isTrimRightBoundary',true,'format','list'));
[X,Y] = get_wavevectors(N_wv,1,struct('isTrimRightBoundary',true,'format','grid'));
X(1,2) - X(1,1);
Y(2,1) - Y(1,1);
abs((X(1,2) - X(1,1)) - (Y(2,1) - Y(1,1)))