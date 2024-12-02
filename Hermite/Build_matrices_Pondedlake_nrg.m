function [I, J, K, X, Y] = Build_matrices_Pondedlake_nrg(Phi, Psi_pi, Psi_w, quadrature_nodes, quadrature_weights, N)

I = diag(sqrt(pi)*2.^(2*(1:N)-1).*factorial(2*(1:N)-1));

J = diag(sqrt(pi)*(-(1:N)).*2.^(2*(1:N)-1).*factorial(2*(1:N)-1));

 K = zeros(N,N);

for k = 0 : N-1           
    for j = 0 : N-1       
        factor = (1/4)*sqrt(pi) * sqrt(2^(2*(k+1)-1)*factorial(2*(k+1)-1))*sqrt(2^(2*(j+1)-1)*factorial(2*(j+1)-1));
        if(j==k)
            K(j+1,k+1) = factor*(2*(k+1)-(1/2));
        elseif(j == k-1)
            K(j+1,k+1) = -(1/2)*factor*sqrt((2*(k+1)-1)*(2*(k+1)-2));
        elseif(j==k+1)
            K(j+1,k+1) = -(1/2)*factor*sqrt((2*(k+1))*(2*(k+1)+1));
        end        
    end
end

X = zeros(N,N);

for k = 0 : N-1    
    
    yk = Psi_pi(:,k+1);
    for j = 0 : N-1       
                        
        yj = Phi(:,j+1);
        X(j+1,k+1) = sum( yk.*yj.*quadrature_weights);        
    end
end

%X = 0.5*(X + X');

Y = zeros(N,N);

for k = 0 : N-1    
    
    yk = Psi_w(:,k+1);
    for j = 0 : N-1       
                        
        yj = Phi(:,j+1);
        Y(j+1,k+1) = sum( exp(-quadrature_nodes.^2/4).*yk.*yj.*quadrature_weights);        
    end
end

Y = 0.5*(Y + Y');


end