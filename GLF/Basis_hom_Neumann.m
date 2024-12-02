function [Psi] = Basis_hom_Neumann(N, x, L)

Psi = zeros(length(x),N);

for p = 0 : N-1    
    Psi(:,p+1) = exp(-0.5*x).*(L(:,p+1) - ((2*p+1)/(2*p+3))*L(:,p+2));    
end

end

