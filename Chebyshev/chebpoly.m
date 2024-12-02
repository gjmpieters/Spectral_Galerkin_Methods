function v = chebpoly(k,x)

%v = (1/2)*( ( x + sqrt(x.^2 - 1) ).^k + ( x - sqrt(x.^2 - 1) ).^k );

v = cos(k*acos(x));
