clear; close all;

fig = figure;

tiledlayout('flow')

for i = 1:4
    for j = 1:4
        nexttile
        ax(i,j) = gca;
        scatter(-1.5 + 3*rand*rand(50,1),-1.5 + 3*rand*rand(50,1))
    end
end