load u_values_jr.txt -ascii

u = u_values;
N = 100;
h = 1/(N+1);

x=0:h:N*h;

plot(x,u)
hold on

%%

load u_values_mh.txt -ascii

u2 = u_values;
N2 = 100;
h2 = 1/(N2+1);

x2=0:h2:N2*h2;

plot(x2,u2,'r--')