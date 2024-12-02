function [I, D, C, E, T, J, G] = Build_matrices_Stuttgart(Phi, Psi, quadrature_nodes, quadrature_weights, n, Der, t)

E = zeros(n,n);

for k = 0 : n-1    
    
    yk = Der*Der*Phi(:,k+1);
    for j = 0 : n-1
                
        yj = Phi(:,j+1);
        E(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end


C = zeros(n,n);

for k = 0 : n-1    
    
    yk = Psi(:,k+1);
    for j = 0 : n-1
                
        yj = Phi(:,j+1);
        C(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end


I = zeros(n,n);

for k = 0 : n-1    
    
    yk = Phi(:,k+1);
    for j = 0 : n-1
                
        yj = Phi(:,j+1);
        I(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end

D = zeros(n,n);

for k = 0 : n-1    
    
    yk = Psi(:,k+1);
    for j = 0 : n-1
                
        yj = Psi(:,j+1);
        D(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end

J = zeros(n,n);

for k = 0 : n-1    
    
    yk = Der*Psi(:,k+1);
    for j = 0 : n-1
                
        yj = Psi(:,j+1);
        J(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end

T = zeros(n,n);

for k = 0 : n-1    
    
    yk = Der*Der*Psi(:,k+1);
    for j = 0 : n-1                             
                
        yj = Psi(:,j+1);
        T(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end

G = zeros(n,n);

for k = 0 : n-1    
    
    yk = Phi(:,k+1);
    g = Groundstate_Stuttgart(quadrature_nodes, t);
    for j = 0 : n-1                             
                
        yj = Psi(:,j+1);
        G(j+1,k+1) = sum( g.*yk.*yj.*quadrature_weights);
    end
end

end

