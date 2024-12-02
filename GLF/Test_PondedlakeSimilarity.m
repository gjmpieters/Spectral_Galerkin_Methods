clear

figure(1)

n = 32;

[quadrature_nodes, quadrature_weights] = Laguerre_Gauss_quadrature(n);
[L] = GenerateLaguerrePolynomials(n, quadrature_nodes);
[Phi] = Basis_hom_Dirichlet(n, quadrature_nodes, L);
[Psi] = Basis_hom_Neumann(n, quadrature_nodes, L);
Der = Derivative(quadrature_nodes, n, L);

[C, D, S, T, B, N, O, E, F, G, H, J, K] = Build_matrices_Pondedlake(Phi, Psi, quadrature_nodes, quadrature_weights, n, Der);

Da = (150.0 - 0.0001)/50000;
a = 0.00001:Da:150.0;

% for i=1:length(a)
%      i;
%      LHS = S + (1/4)*C + (1/16)*F;
%      [E_Cheb, M_Cheb] = solve_evp(inv(C)*LHS,0,'Spectrum');
%  end
%  
% plot(a,R_stab,'y','LineWidth',1.0);
% grid on
% axis([0 3 0 16])
% axis square
 

for i=1:length(a)
    i;
    BB = -T - a(i)^2*D;
    %LHS = -S - a(i)^2*C + (1/2)*E;
    LHS = -S - a(i)^2*C - (1/4)*C - (1/16)*F;
    %RHS = - a(i)^2*N*inv(BB)*B;
    RHS =  a(i)^2*H*inv(BB)*G;
    [E_Cheb, M_Cheb] = solve_evp(inv(LHS)*RHS,0,'Spectral-R-linearised');
    E_Cheb = 1./E_Cheb;
    R_lin(i) = E_Cheb(1);
end

[value idx] = min(R_lin);
a(idx)
R_lin(idx)

plot(a,R_lin,'b','LineWidth',1.0);
grid on
axis([0 3 0 16])
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
