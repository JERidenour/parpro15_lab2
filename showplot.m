load u_values.txt -ascii

u = u_values;
N = 100;
h = 1/(N+1);

x=0:h:N*h;

plot(x,u)