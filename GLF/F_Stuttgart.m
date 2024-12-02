function [f] = F_Stuttgart(quadrature_nodes, tau)

f = 1-(1/2)*exp(-quadrature_nodes).*erfc((quadrature_nodes-tau)./(2*sqrt(tau)))-(1/2)*erfc((quadrature_nodes+tau)./(2*sqrt(tau)));

end
