function f = Dw(k,x)

f = k * sqrt(1 - x.^2) .* sin(k * acos(x)) - 2*x.*chebpoly(k,x);