function [g_der] = Groundstate_Stuttgart(quadrature_nodes, t)

[f] = F_Stuttgart(quadrature_nodes, t);
%[f_der] = F_deriv_Stuttgart(quadrature_nodes, t);

%plot(quadrature_nodes,f)
%hold on
%plot(quadrature_nodes,f_der)
%hold on

fun = @(t)F_deriv_Stuttgart(quadrature_nodes,t);

g = 1 + integral(fun, 0, t, 'ArrayValued',true);
g_der = f - g;

% plot(quadrature_nodes,g)
% hold on
% plot(quadrature_nodes,g_der)
% axis square
% grid on

end
