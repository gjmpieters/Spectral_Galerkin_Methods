clear

N = 32;
alpha = 0;

[quadrature_nodes, quadrature_weights] = Laguerre_Gauss_quadrature(N, alpha);
[Phi] = Basis_hom_Dirichlet(N,quadrature_nodes, 1);
D = Derivative(quadrature_nodes, N);
[I, S, A, M, E] = Build_matrices_Saltlake(Phi, quadrature_nodes, quadrature_weights, N, D);

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

Da = (3.0 - 0.001)/1000;
a = 0.001:Da:3.0;

for i=1:length(a)
    i;
    B = -S - a(i)^2*I;
    LHS = -S - a(i)^2*I - E;
    RHS = - a(i)^2*M*inv(B)*I;
    [E_Cheb, M_Cheb] = solve_evp(inv(LHS)*RHS,0,'Spectral-R-linearised');
    E_Cheb = 1./E_Cheb;
    R_stab(i) = E_Cheb(1);
end

[value idx] = min(R_stab);
a(idx)

plot(a,R_stab,'r--','LineWidth',1.0);
grid on
axis([0 3 0 50])
axis square

hold on

for i=1:length(a)
    i;
    B = -S - a(i)^2*I;
    LHS = -S - a(i)^2*I - 0.25*S;
    RHS = - 0.5*a(i)^2*(M*inv(B)*I + I*inv(B)*M);
    [E_Cheb, M_Cheb] = solve_evp(inv(LHS)*RHS,0,'Spectral-R-energy');
    E_Cheb = 1./E_Cheb;
    R_stab(i) = E_Cheb(1);
end

plot(a,R_stab,'b--','LineWidth',1.0);
grid on
axis([0 3 0 50])
axis square

[value idx] = min(R_stab);
a(idx)