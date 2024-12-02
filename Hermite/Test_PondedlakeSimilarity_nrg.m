clear

figure(1)

n = 5;

%[quadrature_nodes, quadrature_weights] = Hermite_Gauss_quadrature(2*n);
[quadrature_nodes, quadrature_weights] = NewtonCotes(50,100);
[H] = GenerateHermiteFunctions(2*n, quadrature_nodes, 0.125);
[Phi] = Basis_hom_Dirichlet(n, quadrature_nodes, H);

Da = (3.0 - 0.01)/100;
%a = 0.01:Da:3.0;
a = 3.0:-Da:0.01;
%Da = (150.0 - 0.0001)/50000;
%a = 0.00001:Da:150.0;

R_nrg = 17;

for i=1:length(a)
    
    G = Green_function(quadrature_nodes, quadrature_nodes, a(i));    
    [Psi_pi] = Basis_hom_Neumann_pi(n, quadrature_nodes, quadrature_weights, H, G);
    [Psi_w] = Basis_hom_Neumann_w(n, quadrature_nodes, quadrature_weights, H, G);
    [I, J, K, X, Y] = Build_matrices_Pondedlake_nrg(Phi, Psi_pi, Psi_w, quadrature_nodes, quadrature_weights, n);
    
    %conditioner = pinv(X);
    Di = diag(diag(X));
    
    options = optimset('Display','iter'); % show iterations
       
    fun = @(R) det(inv(Di)*(-K - (a(i)^2 + 0.25)*I - R*0.5*a(i)^2*(Y+X)));    
    %fun = @(R) det(inv(Di)*(J - a(i)^2*I - R*a(i)^2*Y));
    [sol,fval] = fzero(fun,R_nrg(end),options)
    
    %test = [2:30];
    %for t = 1 : length(test)
    %    v(t) = fun(test(t));
    %end
    %plot(test, v);
    
    %e = eig(a(i)^2*Y, J - a(i)^2*I);
    
    R_nrg(i) = sol
end

[value idx] = min(R_nrg);
a(idx)
R_nrg(idx)

plot(a,R_nrg,'b','LineWidth',1.0);
grid on
axis([0 3 0 16])
axis square

hold on

% for i=1:length(a)
%     i;
%     BB = -T - a(i)^2*D;
%     LHS = -S - a(i)^2*C - 0.25*C;
%     RHS = - 0.5*a(i)^2*(B'*inv(BB)*O + N*inv(BB)*B);
%     [E_Cheb, M_Cheb] = solve_evp(inv(LHS)*RHS,0,'Spectral-R-energy');
%     E_Cheb = 1./E_Cheb;
%     R_nrg(i) = E_Cheb(1);
% end
% 
% plot(a,R_nrg,'r','LineWidth',1.0);
% grid on
% axis([0 3 0 16])
% axis square
% 
% [value idx] = min(R_nrg);
% a(idx)
% 
% plot(a, (3*pi/32)*(a.^2 + 1).*(8*a.^2+1)./a.^2, 'g','LineWidth',1.0);
% 
% hold on

xlabel('$$b$$','Interpreter','latex');
ylabel('$$R^{*}$$','Interpreter','latex');

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
