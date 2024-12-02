% function [Phi] = Basis_hom_Dirichlet(N, x, beta)
% 
% Phi = zeros(length(x),N);
% 
% parfor p = 0 : N-1    
%     Phi(:,p+1) = exp(-0.5*beta*x).*(laguerreL(p,beta*x) - laguerreL(p+1,beta*x));    
% end
% 
% end
% 

function [Phi] = Basis_hom_Dirichlet(N, x, L)

Phi = zeros(length(x),N);

for p = 0 : N-1
    Phi(:,p+1) = exp(-0.5*x).*(L(:,p+1) - L(:,p+2));    
end

end

