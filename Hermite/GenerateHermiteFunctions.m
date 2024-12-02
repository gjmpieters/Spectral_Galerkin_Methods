function [H] = GenerateHermiteFunctions(N,x,coef)

H = zeros(length(x),N);

H(:,1) = 1*exp(-x.^2*coef);
H(:,2) = x.*exp(-x.^2*coef);

for n = 1 : N-1  
  H(:,n+2) = x.*H(:,n+1) - 2*n*H(:,n);  
end




