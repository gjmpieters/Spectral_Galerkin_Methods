function [quadrature_nodes, quadrature_weights] = Laguerre_Gauss_quadrature(N)

% Build matrix (7.38) as in Shen
a = 2 * (0:N) + 1;
b = (1:N).*((1:N));

A = diag(a) + diag(-sqrt(b),1) + diag(-sqrt(b),-1);

quadrature_nodes = eig(A);

quadrature_weights = (1 / ((N+1)^2)) * quadrature_nodes ./ (Laguerre_function(N,0,quadrature_nodes).*Laguerre_function(N,0,quadrature_nodes));

%else
%    quadrature_weights = (gamma(N + alpha + 1) / ((N+alpha+1)*factorial(N+1))) * quadrature_nodes ./ (Laguerre_function(N,alpha,quadrature_nodes).*Laguerre_function(N,alpha,quadrature_nodes));
%end


