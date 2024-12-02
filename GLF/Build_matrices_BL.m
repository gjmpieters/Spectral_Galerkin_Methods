function [mI, mE, mF, mG, mH, mJ_I, mJ_E, mJ_F, mJ_G, mJ_H] = Build_matrices_BL(Phi, quadrature_nodes, quadrature_weights, D, N, s, no, nw, M, lambda)

f_prime = F_prime(min(1,max(eps,Phi*s)), no, nw, M);
f_primeprime = F_primeprime(min(1,max(eps,Phi*s)), no, nw, M);
%f_prime = F_prime(Phi*s, no, nw, M);
%f_primeprime = F_primeprime(Phi*s, no, nw, M);
g = G(min(1,max(eps,Phi*s)), no, nw, M, lambda);
g_prime = G_prime(min(1,max(eps,Phi*s)), no, nw, M, lambda);

mJ_I = zeros(N,N);

for k = 0 : N-1    
    
    yk = Phi(:,k+1);
    for j = 0 : N-1
                
        yj = Phi(:,j+1);
        mJ_I(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end

mI = mJ_I*s;


mJ_E = zeros(N,N);

for k = 0 : N-1    
    
    yk = D*Phi(:,k+1);
    for j = 0 : N-1
                
        yj = Phi(:,j+1);
        mJ_E(j+1,k+1) = sum( quadrature_nodes.*yk.*yj.*quadrature_weights);
    end
end

mE = mJ_E*s;

mFtmp = zeros(N,N);

for k = 0 : N-1    
    
    yk = D*Phi(:,k+1);
    for j = 0 : N-1       
                        
        yj = Phi(:,j+1);
        mFtmp(j+1,k+1) = sum( f_prime.*yk.*yj.*quadrature_weights);
    end
end

mF = mFtmp*s;

mGtmp = zeros(N,N);

for k = 0 : N-1    
    
    yk = D*Phi(:,k+1);
    for j = 0 : N-1                             
                
        yj = D*Phi(:,j+1);
        mGtmp(j+1,k+1) = sum( quadrature_nodes.*g.*yk.*yj.*quadrature_weights);                
    end
end

mG = mGtmp*s;

mJ_H = zeros(N,N);

for k = 0 : N-1    
    
    yk = D*D*Phi(:,k+1);
    for j = 0 : N-1
                
        yj = Phi(:,j+1);
        mJ_H(j+1,k+1) = sum( yk.*yj.*quadrature_weights);
    end
end

mH = mJ_H*s;


mJ_F = zeros(N,N);
mJ_Ftmp = zeros(1,N);

for i = 0 : N-1    
    
    yi = Phi(:,i+1);
    dyi = D*Phi(:,i+1);
    for j = 0 : N-1    
    
        yj = Phi(:,j+1);
        for k = 0 : N-1       
        
            yk = D*Phi(:,k+1);
            mJ_Ftmp(1,k+1) = sum( f_primeprime.*yi.*yk.*yj.*quadrature_weights);
        end
                
        mJ_F(j+1, i+1) = mJ_Ftmp*s + sum( f_prime.*dyi.*yj.*quadrature_weights);
    end    
end

mJ_G = zeros(N,N);
mJ_Gtmp = zeros(1,N);

for i = 0 : N-1    
    
    yi = Phi(:,i+1);
    dyi = D*Phi(:,i+1);
    for j = 0 : N-1    
    
        yj = D*Phi(:,j+1);
        for k = 0 : N-1       
        
            yk = D*Phi(:,k+1);
            mJ_Gtmp(1,k+1) = sum( quadrature_nodes.*g_prime.*yi.*yk.*yj.*quadrature_weights);            
        end
                
        mJ_G(j+1, i+1) = mJ_Gtmp*s + sum( quadrature_nodes.*g.*dyi.*yj.*quadrature_weights);        
    end    
end

end

