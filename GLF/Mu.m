function [m] = Mu(quadrature_nodes)

m = exp(-quadrature_nodes.^2/8);

end
