function [f_der] = F_deriv_Stuttgart(quadrature_nodes, tau)

if(tau == 0)
    f_der= 0;
else
    f_der = (1/2)*exp(-quadrature_nodes).*erfc((1/2)*(quadrature_nodes-tau)./sqrt(tau))+(1/2)*exp(-quadrature_nodes).*exp(-(1/4)*(quadrature_nodes-tau).^2./tau)./(sqrt(pi)*sqrt(tau))+(1/2)*exp(-(1/4)*(quadrature_nodes+tau).^2./tau)./(sqrt(pi)*sqrt(tau));
end

end
