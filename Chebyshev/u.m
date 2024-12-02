function f = u(k,x)

c = (4*k+4)/(-2*k^2-6*k-5);
d = (2*k^2+2*k+1)/(-2*k^2-6*k-5);

f = chebpoly(k,x) + c * chebpoly(k+1,x) + d * chebpoly(k+2,x);