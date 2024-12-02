function [y] = Laguerre_function(N, alpha, x)

y = exp(-x./2).*laguerreL(N,alpha,x);

end

