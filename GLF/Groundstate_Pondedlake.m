function [g] = Groundstate_Pondedlake(quadrature_nodes)

g = -exp(-0.25*quadrature_nodes.^2);

end
