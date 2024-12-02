function f = Dv(k,x)

f = (k*sin(k*acos(x)) - (k+2)*sin((k+2)*acos(x))) ./ sqrt(1 - x.^2);