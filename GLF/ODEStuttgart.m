function dy = ODEStuttgart(t, y, quadrature_nodes)

[f, fder] = F_Stuttgart(quadrature_nodes, t);

dy = fder;

end