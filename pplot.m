load u_values_jr.txt -ascii

u = u_values_jr;
N = 100;
h = 1/(N+1);

x=0:h:N*h;

plot(x,u)
hold on

%%

load u_values_mh.txt -ascii

u2 = u_values_mh;
N2 = 100;
h2 = 1/(N2+1);

x2=0:h2:N2*h2;

plot(x2,u2,'r--')

%%

load u_values_par.txt -ascii

u3 = u_values_par;
N3 = length(u3)-1;
h3 = 1/(N3);

x3=0:h3:N3*h3;

plot(x3,u3,'r--')