clear

figure(1)

n = 32;
t = 0.3;
%tau = 5;

[quadrature_nodes, quadrature_weights] = Laguerre_Gauss_quadrature(n);
[L] = GenerateLaguerrePolynomials(n, quadrature_nodes);
[Phi] = Basis_hom_Dirichlet(n, quadrature_nodes, L);
[Psi] = Basis_hom_Robin(n, quadrature_nodes, L);
Der = Derivative(quadrature_nodes, n, L);

aMax = 3;
Da = (aMax - 0.001)/100;
a = 0.001:Da:aMax;

%[I, D, C, E, T, J, G] = Build_matrices_Stuttgart(Phi, Psi, quadrature_nodes, quadrature_weights, n, Der, tau);    

x0 = [0.00078, 0.01668]; %t=100
%x0 = [0.0078, 0.1668]; %t=10
%x0 = [0.01, 2]; %t=1,2
%x0 = [0.01, 7]; %t=0.5

for i=1:length(a)
    i
    
    options = optimset('Display','iter', 'TolFun', 1e-6, 'TolX', 1e-6); % show iterations
    fun = @(x) Stuttgart_goalseek(n, t, a(i), x);
    res = fminsearch(fun, x0, options);
    
    ahat(i) = res(1);
    tau(i) = res(2);
    x0 = [ahat(i), tau(i)];
    
    test = Stuttgart_goalseek(n,t,a(i),x0);
    
    %tau_interval = [max(0.001,tau(i)-2.5) tau(i)+4.5];
    
    [I, D, C, E, T, J, G] = Build_matrices_Stuttgart(Phi, Psi, quadrature_nodes, quadrature_weights, n, Der, tau(i));    
    
    BB = E - ahat(i)^2*I;
    LHS = T + J - ahat(i)^2*D;
    RHS =  -ahat(i)^2*G*inv(BB)*C;
    [E_Cheb, M_Cheb] = solve_evp(inv(LHS)*RHS,0,'Spectral-R-linearised');
    E_Cheb = 1./E_Cheb;
    R_lin(i) = E_Cheb(1,1);
    
    R_E(i) = 1/R_lin(i);
    
    %     z = 0:0.01:10;
    %     [L_tmp] = GenerateLaguerrePolynomials(n, z');
    %     [Phi_tmp] = Basis_hom_Dirichlet(n, z', L_tmp);
    %     [Psi_tmp] = Basis_hom_Robin(n, z', L_tmp);
    %     sol_x = BuildSolution(M_Cheb(:,1), Psi_tmp);
    %     sol_w = BuildSolution(-a(i)^2*inv(BB)*C*M_Cheb(:,1), Phi_tmp);
    %     plot(z,sol_x);
    %     hold on
    %     plot(z,sol_w);
    %     clf(1)
    %     Groundstate_Stuttgart(z,0.1);
    %     Groundstate_Stuttgart(z,tau);
    
    
end

% [value idx] = min(R_lin);
% a(idx)
% R_lin(idx)

%figure(2)
plot(a,R_E,'b','LineWidth',1.0);
grid on
axis([0 3 0 10])
axis square

hold on


for i=1:length(a)
    i;
    BB = -T - a(i)^2*D;
    LHS = -S - a(i)^2*C - 0.25*C;
    RHS = - 0.5*a(i)^2*(B'*inv(BB)*O + N*inv(BB)*B);
    [E_Cheb, M_Cheb] = solve_evp(inv(LHS)*RHS,0,'Spectral-R-energy');
    E_Cheb = 1./E_Cheb;
    R_nrg(i) = E_Cheb(1);
end

plot(a,R_nrg,'r','LineWidth',1.0);
grid on
axis([0 3 0 16])
axis square

[value idx] = min(R_nrg);
a(idx)

plot(a, (3*pi/32)*(a.^2 + 1).*(8*a.^2+1)./a.^2, 'g','LineWidth',1.0);

hold on

xlabel('$$\hat{k}$$','Interpreter','latex');
ylabel('$$\hat{R}$$','Interpreter','latex');

line([0,3],[0,12.775],'Color','black')
line([0,3],[0,13.65],'Color','black')

figure(2)
spl = @(b) spline(a,R_nrg,b);
Rmax = 7.56;
Rrange = [Rmax 15 20];
aMinrange = [0.1 0.1 0.1];
numpoints = 200;
for i=1:length(Rrange)
    Rs = Rrange(i);
        
    deltatesta = (Rs/Rmax - aMinrange(i)) / numpoints;
    testa = aMinrange(i):deltatesta:Rs/Rmax;
    alpha = Rs ./ (testa * sqrt(pi));
    
    for j=1:length(alpha)
        %for i=1:1
        alph = alpha(j)
        %plot(a,alph*a);
        
        fun = @(b) spl(b) - alph*b;
                       
        bstarL(j) = fzero(fun,[a(1) 1.25]);       
        bstarR(j) = fzero(fun,[1.25 a(end)]);        
    end        
   
    loglog((bstarL ./ testa(1:end)).^2, testa(1:end),'r--')
    hold on
    loglog((bstarR ./ testa(1:end)).^2, testa(1:end),'b--')    
end

axis([0.01 1000 0.5 4])
axis square
grid on
xlabel('$$t$$','Interpreter','latex');
ylabel('$$a$$','Interpreter','latex');    

%figure(3)
spl = @(b) spline(a,R_lin,b);
Rmax = 8.08;
Rrange = [7.56 15 20];
aMinrange = [0.1 0.1 0.1];
numpoints = 200;
for i=1:length(Rrange)
    Rs = Rrange(i);
        
    deltatesta = (Rs/Rmax - aMinrange(i)) / numpoints;
    testa = aMinrange(i):deltatesta:Rs/Rmax;
    alpha = Rs ./ (testa * sqrt(pi));
    
    for j=1:length(alpha)
        %for i=1:1
        alph = alpha(j)
        %plot(a,alph*a);
        
        fun = @(b) spl(b) - alph*b;
                       
        bstarL(j) = fzero(fun,[a(1) 1.5]);       
        bstarR(j) = fzero(fun,[1.5 a(end)]);        
    end        
   
    loglog((bstarL ./ testa(1:end)).^2, testa(1:end),'r:')
    hold on
    loglog((bstarR ./ testa(1:end)).^2, testa(1:end),'b:')    
end

axis([0.01 1000 0.5 4])
axis square
grid on
xlabel('$$t$$','Interpreter','latex');
ylabel('$$a$$','Interpreter','latex');    
