function [Psi] = Basis_hom_Neumann_pi(N, x, w, H, G)

Psi = zeros(length(x),N);

for p = 0 : N-1    
    
    %for i=1:length(x)                     
    %    sum = 0;
    %    for j=1:length(x)
    %        sum = sum + G(i,j) * H(j, 2*(p+1)) * exp(-x(j)^2/8) * w(j);
    %    end
    %  
    %    Psi(i,p+1) = sum;
    %end 
    
    Psi(:,p+1) = sum(G(:,:)' .* H(:, 2*(p+1)) .* exp(-x.^2/4) .* w);    
end

end

