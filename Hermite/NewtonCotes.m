function [quadrature_nodes, quadrature_weights] = NewtonCotes(N, L)

h = 1/N;
quadrature_nodes = [0:h:L]';
quadrature_weights = h * ones(length(quadrature_nodes),1);

