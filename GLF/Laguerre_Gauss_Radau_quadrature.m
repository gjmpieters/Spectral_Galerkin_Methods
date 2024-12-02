function [quadrature_nodes, quadrature_weights] = Laguerre_Gauss_Radau_quadrature(N, alpha)

% Build matrix (7.38) as in Shen. Note interior nodes are the zeros of the
% derivative
alpha = alpha + 1;
N = N - 1;

a = 2 * (0:N) + alpha + 1;
b = (1:N).*((1:N) + alpha);

A = diag(a) + diag(-sqrt(b),1) + diag(-sqrt(b),-1);

quadrature_nodes = eig(A);
quadrature_weights = (gamma(N + alpha + 1) / ((N+alpha+1)*factorial(N+1))) * quadrature_nodes ./ (Laguerre_function(N,alpha,quadrature_nodes).*Laguerre_function(N,alpha,quadrature_nodes));

