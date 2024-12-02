mfunction f = D2u(k,x)

c = (4*k+4)/(-2*k^2-6*k-5);
d = (2*k^2+2*k+1)/(-2*k^2-6*k-5);

f = ( -k^2*chebpoly(k,x) + x.*Du(k,x) + c*(-(k+1)^2*chebpoly(k+1,x)) + d*(-(k+2)^2*chebpoly(k+2,x)) ) ./ (1 - x.^2);