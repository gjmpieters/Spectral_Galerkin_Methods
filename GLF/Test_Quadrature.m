clear

N = 5;
alpha = 0;
[quadrature_nodes, quadrature_weights] = Laguerre_Gauss_quadrature(N,alpha);
%[quadrature_nodesR, quadrature_weightsR] = Laguerre_Gauss_Radau_quadrature(N,alpha);
%quadrature_nodesR = [0 quadrature_nodesR']';

D = Derivative(N);

uq = Laguerre_function(N,alpha,quadrature_nodes);
%uR = Laguerre_function(N+1,alpha,quadrature_nodesR);

plot(quadrature_nodes,D*uq, 'r*')
hold on
plot(quadrature_nodes,D*D*uq, 'bo')
hold on
%plot(quadrature_nodesR,0, 'o')
%hold on

x = 0:0.001:50;

u = Laguerre_function(N,alpha,x);
uprime = -Laguerre_function(N-1,alpha+1,x) - 0.5 * Laguerre_function(N,alpha,x);
uprimeprime = Laguerre_function(N-2,alpha+2,x) + Laguerre_function(N-1,alpha+1,x) + 0.25 * Laguerre_function(N,alpha,x);

plot(x,u,'k')
hold on
plot(x,uprime,'r')
plot(x,uprimeprime,'b')

axis([0 60 -3 3])
grid on


sum = 0;
for i = 1 : length(x)-1
    delta = x(i+1) - x(i);
       
    sum = sum + 0.5*(u(i)*uprime(i) + u(i+1)*uprime(i+1))*delta;
end
sum