function [P] = Build_initialcondition_Stuttgart(Psi, quadrature_nodes, quadrature_weights, n, R)
         
P = zeros(n,1);

IC = 0.001*exp(-2*R*quadrature_nodes);

for j = 0 : n-1       
    yj = Psi(:,j+1);
    P(j+1) = sum( IC.*yj.*quadrature_weights);
end

end

