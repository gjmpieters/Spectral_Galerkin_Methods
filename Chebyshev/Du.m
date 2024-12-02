function f = Du(k,x)

c = (4*k+4)/(-2*k^2-6*k-5);
d = (2*k^2+2*k+1)/(-2*k^2-6*k-5);

f = (k*sin(k*acos(x)) + c * (k+1)*sin((k+1)*acos(x)) + d * (k+2)*sin((k+2)*acos(x))) ./ sqrt(1 - x.^2);