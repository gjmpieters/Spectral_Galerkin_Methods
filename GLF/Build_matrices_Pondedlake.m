function [C, D, S, T, B, N, O, E, F, G, H, J, K] = Build_matrices_Pondedlake(Phi, Psi, quadrature_nodes, quadrature_weights, n, Der)
         
E = zeros(n,n);

for k = 0 : n-1    
    
    yk = Der*Phi(:,k+1);
    for j = 0 : n-1
                
        yj = Phi(:,j+1);
        xi = quadrature_nodes;
        E(j+1,k+1) = sum( xi.*yk.*yj.*quadrature_weights);
    end
end

F = zeros(n,n);

for k = 0 : n-1    
    
    yk = Phi(:,k+1);
    for j = 0 : n-1
                
        yj = Phi(:,j+1);
        xi_squared = quadrature_nodes.*quadrature_nodes;
        F(j+1,k+1) = sum( xi_squared.*yk.*yj.*quadrature_weights);
    end
end

K = zeros(n,n);

for k = 0 : n-1    
    
    yk = Psi(:,k+1);
    for j = 0 : n-1
                
        yj = Psi(:,j+1);
        xi_squared = quadrature_nodes.*quadrature_nodes;
        K(j+1,k+1) = sum( xi_squared.*yk.*yj.*quadrature_weights);
    end
end

G = zeros(n,n);

for k = 0 : n-1    
    
    yk = Phi(:,k+1);
    for j = 0 : n-1
                
        yj = Psi(:,j+1);
        m = Mu(quadrature_nodes);        
        G(j+1,k+1) = sum( m.*yk.*yj.*quadrature_weights);
    end
end

H = zeros(n,n);

for k = 0 : n-1    
    
    yk = Psi(:,k+1);
    for j = 0 : n-1       
                        
        yj = Phi(:,j+1);
        m = Mu(quadrature_nodes);       
        H(j+1,k+1) = sum( m.*yk.*yj.*quadrature_weights);
    end
end

C = zeros(n,n);

for k = 0 : n-1    
    
    yk = Phi(:,k+1);
    for j = 0 : n-1
                
        yj = Phi(:,j+1);
        C(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
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
        xi = quadrature_nodes;
        J(j+1,k+1) = sum( xi.*yk.*yj.*quadrature_weights);
    end
end

N = zeros(n,n);

for k = 0 : n-1    
    
    yk = Psi(:,k+1);
    for j = 0 : n-1       
                        
        yj = Phi(:,j+1);
        g = Groundstate_Pondedlake(quadrature_nodes);
        N(j+1,k+1) = sum( g.*yk.*yj.*quadrature_weights);
    end
end

O = zeros(n,n);

for k = 0 : n-1    
    
    yk = Phi(:,k+1);
    for j = 0 : n-1       
                        
        yj = Psi(:,j+1);
        g = Groundstate_Pondedlake(quadrature_nodes);
        O(j+1,k+1) = sum( g.*yk.*yj.*quadrature_weights);
    end
end

S = zeros(n,n);

for k = 0 : n-1    
    
    yk = Der*Phi(:,k+1);
    for j = 0 : n-1                             
                
        yj = Der*Phi(:,j+1);
        S(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end

T = zeros(n,n);

for k = 0 : n-1    
    
    yk = Der*Psi(:,k+1);
    for j = 0 : n-1                             
                
        yj = Der*Psi(:,j+1);
        T(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end

B = zeros(n,n);

for k = 0 : n-1    
    
    yk = Phi(:,k+1);
    for j = 0 : n-1       
                        
        yj = Psi(:,j+1);
        B(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end


end

