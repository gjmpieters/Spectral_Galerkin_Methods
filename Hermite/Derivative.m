function [D] = Derivative(quadrature_nodes)

N = length(quadrature_nodes);
delta = quadrature_nodes(2) - quadrature_nodes(1);
D = diag(ones(N-1,1)/(2*delta),+1) + diag(-ones(N-1,1)/(2*delta),-1);
D(1,2) = 1/delta;
D(1,1) = -1/delta;

D(N,N) = 1/delta;
D(N,N-1) = -1/delta;

