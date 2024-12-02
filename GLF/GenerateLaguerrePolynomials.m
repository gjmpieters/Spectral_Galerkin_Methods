function [L] = GenerateLaguerrePolynomials(N,x)

L = zeros(length(x),N);

L(:,1) = 1;
L(:,2) = (1-x);

for n = 1 : N-1
  %test = laguerreL(n+1,x);
  %plot(test)
  %hold on
  
  L(:,n+2) = ((2*n+1-x).*L(:,n+1) - n*L(:,n))/(n+1);
  %plot(L(:,n+2))
end




