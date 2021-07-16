clear; close all;
f = @(n) n.*(ceil(n/2));

n_c = 1:0.01:10;

plot(n_c,f(n_c))

n_i = 1:10;

hold on
plot(n_i,f(n_i),'k*')

legend('continuous','integer')
ylabel('num\_sample\_points')
xlabel('N\_wv(1)')