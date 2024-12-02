function [I, J, X] = Build_matrices_Pondedlake(Phi, Psi, quadrature_nodes, quadrature_weights, N)

I = diag(sqrt(pi)*2.^(2*(1:N)-1).*factorial(2*(1:N)-1));

J = diag(sqrt(pi)*(-(1:N)).*2.^(2*(1:N)-1).*factorial(2*(1:N)-1));

% I = zeros(N,N);
% 
% for k = 0 : N-1    
%     
%     yk = Phi(:,k+1);
%     for j = 0 : N-1
%                 
%         yj = Phi(:,j+1);
%         I(j+1,k+1) = sum( yk.*yj.*quadrature_weights);        
%     end
% end

% J = zeros(N,N);
% 
% for k = 0 : N-1    
%     
%     yk = Phi(:,k+1);
%     for j = 0 : N-1
%                 
%         yj = Phi(:,j+1);
%         J(j+1,k+1) = (-(j+1))*sum( yk.*yj.*quadrature_weights);
%     end
% end

X = zeros(N,N);

for k = 0 : N-1    
    
    yk = Psi(:,k+1);
    for j = 0 : N-1       
                        
        yj = Phi(:,j+1);
        X(j+1,k+1) = sum( exp(-quadrature_nodes.^2/8).*yk.*yj.*quadrature_weights);        
    end
end

X = 0.5*(X + X');


end