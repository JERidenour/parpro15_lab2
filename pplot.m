%% 

clear
clf



%%
% load u_values_jr.txt -ascii
% 
% u = u_values_jr;
% N = length(u_values_jr)
% h = 1/(N+1);
% 
% x=0:h:(N-1)*h;
% 
% plot(x',u,'b:')
% hold on
% axis tight

%%

load u_values_mh.txt -ascii

u2 = u_values_mh;
N2 = length(u_values_mh)
h2 = 1/(N2-1)

x2=0:h2:1;  % x: [0 , 1]

plot(x2,u2,'m');
hold on
axis tight

%%

load u_values_gs.txt -ascii

ugs = u_values_gs;
Ngs = length(ugs)
hgs = 1/(Ngs-1)

xgs=0:hgs:1;  % x: [0 , 1]

plot(xgs,ugs,'k--');
hold on
axis tight

%%

load u_values_par.txt -ascii

u3 = u_values_par;
N3 = length(u3)-1;
h3 = 1/(N3);


x3=0:h3:N3*h3;

plot(x3,u3,'k--')
hold on
axis tight

%%

load u_values_par_N1000_i1000.txt -ascii

u = u_values_par_N1000_i1000;
N = length(u)-1;
h = 1/(N);


x=0:h:N*h;

plot(x,u,'--')
hold on


load u_values_par_N1000_i10000.txt -ascii

u = u_values_par_N1000_i10000;
N = length(u)-1;
h = 1/(N);


x=0:h:N*h;

plot(x,u,'--')
hold on

load u_values_par_N1000_i100000.txt -ascii

u = u_values_par_N1000_i100000;
N = length(u)-1;
h = 1/(N);


x=0:h:N*h;

plot(x,u,'--')
hold on

load u_values_par_N1000_i1000000.txt -ascii

u = u_values_par_N1000_i1000000;
N = length(u)-1;
h = 1/(N);


x=0:h:N*h;

plot(x,u,'k--')
hold on

load u_values_par_N1000_i10000000.txt -ascii

u = u_values_par_N1000_i10000000;
N = length(u)-1;
h = 1/(N);


x=0:h:N*h;

plot(x,u,'m')
hold on

axis tight
title('Convergence N=1000');
legend('10e3','10e4','10e5','10e6','10e7','location','southwest');

xlabel('x')
ylabel('u(x)')
