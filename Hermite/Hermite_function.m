function [y] = Hermite_function(N, x)

y = exp(-x.^2/8).*hermiteH(N,x);

end

