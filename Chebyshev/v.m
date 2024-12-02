function f = V(k,x)

f = chebpoly(k,x) - chebpoly(k+2,x);