x = 0:0.01:0.3;
y = 100000000*(exp(5*x)-1);
z = -3/4*log(1-4/3*x);

plot(y,x);
figure
plot(y,z);