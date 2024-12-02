clear

n = 32;

[quadrature_nodes, quadrature_weights] = Laguerre_Gauss_quadrature(n);
[L] = GenerateLaguerrePolynomials(n, quadrature_nodes);
[Phi] = Basis_hom_Dirichlet(n, quadrature_nodes, L);
[Psi] = Basis_hom_Neumann(n, quadrature_nodes, L);
Der = Derivative(quadrature_nodes, n, L);
[C, D, S, T, B, N, O, E] = Build_matrices_Pondedlake(Phi, Psi, quadrature_nodes, quadrature_weights, n, Der);

%x = 0:0.1:350;
%f = BuildSolution(x, inv(I)*F);
%plot(x,f,x,2*exp(-x))

% ucoef = inv(S + I)*F;
% 
% x = 0:0.1:350;
% u = BuildSolution(x, ucoef);
% 
% plot(x,u,x,x.*exp(-x),'*')
% axis([0 50 -.5 .5])

Da = (3.0 - 0.0001)/2000;
a = 0.00001:Da:3.0;

for i=1:length(a)
    i;
    BB = -T - a(i)^2*D;   
    LHS = -S - a(i)^2*C;
    RHS = - a(i)^2*N*inv(BB)*B;
    [E_Cheb, M_Cheb] = solve_evp(inv(LHS)*RHS,0,'Spectral-R-linearised');
    E_Cheb = 1./E_Cheb;
    R_stab(i) = E_Cheb(1);
end

[value idx] = min(R_stab);
a(idx)

%plot(a,R_stab,'b--','LineWidth',1.0);
grid on
axis([0 3 0 16])
axis square

hold on

for i=1:length(a)
    i;
    BB = -T - a(i)^2*D;
    LHS = -S - a(i)^2*C;
    RHS = - 0.5*a(i)^2*(B'*inv(BB)*O + N*inv(BB)*B);
    [E_Cheb, M_Cheb] = solve_evp(inv(LHS)*RHS,0,'Spectral-R-energy');
    E_Cheb = 1./E_Cheb;
    R_stab(i) = E_Cheb(1);
end

plot(a,R_stab,'r--','LineWidth',1.0);
grid on
axis([0 3 0 16])
axis square

[value idx] = min(R_stab);
a(idx)

xlabel('$$\hat{k}$$','Interpreter','latex');
ylabel('$$\hat{R}$$','Interpreter','latex');