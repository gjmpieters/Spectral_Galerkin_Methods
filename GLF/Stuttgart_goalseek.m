function F = Stuttgart_goalseek(n, t, a, x)

b = x(1);
tau = x(2);

[quadrature_nodes, quadrature_weights] = Laguerre_Gauss_quadrature(n);
[L] = GenerateLaguerrePolynomials(n, quadrature_nodes);
[Phi] = Basis_hom_Dirichlet(n, quadrature_nodes, L);
[Psi] = Basis_hom_Robin(n, quadrature_nodes, L);
Der = Derivative(quadrature_nodes, n, L);

[I, D, C, E, T, J, G] = Build_matrices_Stuttgart(Phi, Psi, quadrature_nodes, quadrature_weights, n, Der, tau);

BB = E - b^2*I;
LHS = T + J - b^2*D;
RHS =  -b^2*G*inv(BB)*C;
[E_Cheb, M_Cheb] = solve_evp(inv(LHS)*RHS,0,'Spectral-R-linearised');
E_Cheb = 1./E_Cheb;
R_lin = E_Cheb(1,1);


R_lin_desired = sqrt(t/tau);

F = (R_lin - R_lin_desired)^2 + (b - R_lin_desired * a)^2;

end
