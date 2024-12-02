function [quadrature_nodes, quadrature_weights] = Hermite_Gauss_quadrature(N)

% Build matrix (7.38) as in Shen
a = 0 * (0:N);
b = 0.5*(1:N);

A = diag(a) + diag(sqrt(b),1) + diag(sqrt(b),-1);

quadrature_nodes = eig(A);

quadrature_nodes = quadrature_nodes(ceil(end/2):end);

quadrature_weights = (sqrt(pi) * 2^N * factorial(N)) ./ (Hermite_function(N,quadrature_nodes).*Hermite_function(N,quadrature_nodes)) / (N+1);