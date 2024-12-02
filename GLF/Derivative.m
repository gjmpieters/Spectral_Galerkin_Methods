% function [D] = Derivative(quadrature_nodes, N)
% 
% D = zeros(N,N);
% 
% for k = 0 : N
%     for j = 0 : N
%         
%         if(k==j)
%             D(k+1,j+1) = -1 / (2*quadrature_nodes(k+1));
%         else
%             xj = quadrature_nodes(j+1);
%             xk = quadrature_nodes(k+1);
%             D(k+1,j+1) = ((xj * Laguerre_function(N,0,xk)) / (xk * Laguerre_function(N,0,xj))) / (xk - xj);            
%         end
%     end
% end

function [D] = Derivative(quadrature_nodes, N, L)

D = zeros(N,N);

for k = 0 : N
    for j = 0 : N
        
        if(k==j)
            D(k+1,j+1) = -1 / (2*quadrature_nodes(k+1));
        else
            xj = quadrature_nodes(j+1);
            xk = quadrature_nodes(k+1);
            
            wj = exp(-xj/2).*L(j+1,N+1);
            wk = exp(-xk/2).*L(k+1,N+1);
            
            D(k+1,j+1) = ((xj * wk) / (xk * wj)) / (xk - xj);            
            
            %test = ((xj * Laguerre_function(N,0,xk)) / (xk * Laguerre_function(N,0,xj))) / (xk - xj);            
        end
    end
end

