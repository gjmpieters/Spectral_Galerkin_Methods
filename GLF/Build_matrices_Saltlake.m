function [I, S, A, M, E] = Build_matrices_Saltlake(Phi, quadrature_nodes, quadrature_weights, N, D)

E = zeros(N,N);

for k = 0 : N-1    
    
    yk = D*Phi(:,k+1);
    for j = 0 : N-1
                
        yj = Phi(:,j+1);
        E(j+1,k+1) = sum( -0.5*quadrature_nodes.*yk.*yj.*quadrature_weights);
    end
end

I = zeros(N,N);

for k = 0 : N-1    
    
    yk = Phi(:,k+1);
    for j = 0 : N-1
                
        yj = Phi(:,j+1);
        I(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end

M = zeros(N,N);

for k = 0 : N-1    
    
    yk = Phi(:,k+1);
    for j = 0 : N-1       
                        
        yj = Phi(:,j+1);
        g = Groundstate_Saltlake(quadrature_nodes);
        M(j+1,k+1) = sum( g.*yk.*yj.*quadrature_weights);
    end
end

S = zeros(N,N);

for k = 0 : N-1    
    
    yk = D*Phi(:,k+1);
    for j = 0 : N-1                             
                
        yj = D*Phi(:,j+1);
        S(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end

A = zeros(N,N);

for k = 0 : N-1    
    
    yk = D*Phi(:,k+1);
    for j = 0 : N-1       
                        
        yj = Phi(:,j+1);
        A(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end


end

