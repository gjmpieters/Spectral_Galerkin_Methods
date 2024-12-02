function [g] = Groundstate_Saltlake(quadrature_nodes)

%g = -exp(-quadrature_nodes);
g = -exp(-quadrature_nodes.^2/4);

end
