function [Phi] = Basis_hom_Neumann(N, x, H)

Phi = zeros(length(x),N);

for p = 0 : N-1
    Phi(:,p+1) = H(:,2*(p+1));    
end

end

