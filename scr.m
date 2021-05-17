clear; close all;

N = 100;
N_col = 3;
y = rand(N,N_col);
y(:,2) = y(:,2) + .5;
y(:,3) = y(:,3) + 1;

cats = categorical({'cat A','cat B','cat C'});

cats = reshape(repmat(cats,size(y,1),1),[],1);

y = reshape(y,[],1);

figure
b = boxchart(cats,y);
ax = gca();

