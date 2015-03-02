%%
clear
load u_values.txt -ascii

u = u_values;
N = 100;
h = 1/(N+1);

x=0:h:N*h;

plot(x,u)

hold on

%%


N = 1000;
h = 1/(N+1);
tol= 1e-6;
max_diff = 1;

x=0:h:N*h;
f = -x.^(-1.5);
r = x.^3;

%initial guess
u = ones(N+1,1);
u(1)=0;
u(N+1)=0;

while max_diff > tol
    %display(max_diff);
    %pause(1)
   prev_u=u;
   for n = 2:N
        u(n) = (u(n-1)+u(n+1)-h*h*f(n))/(2-h*h*r(n));
   end
   max_diff = max(abs(u-prev_u)); 
end

plot(x,u)