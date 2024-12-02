function [s_sol,w_sol,Dw_sol,D2w_sol,der] = build_Cheb_sol(x,K,H,M_s,M_w)

s_sol = 0;
for i = 0:K
  s_sol = s_sol + M_s(i+1)*w(i,x);
end

w_sol = 0;
for i = 0:K
  w_sol = w_sol + M_w(i+1)*u(i,x);
end

Dw_sol = 0;
for i = 0:K
  Dw_sol = Dw_sol + M_w(i+1)*Du(i,x);
end

der = 0;
D2w_sol = 0;
for i = 0:K
  c = (4*i+4)/(-2*i^2-6*i-5);
  d = (2*i^2+2*i+1)/(-2*i^2-6*i-5);

  der = der + (4/H^2)* M_w(i+1)*((-1)^(i)*((i)^4 + (i)^2)/3 + c * ((-1)^(i+1)*((i+1)^4 + (i+1)^2)/3) + d * ((-1)^(i+2)*((i+2)^4 + (i+2)^2)/3));
  
  D2w_sol = D2w_sol + (4/H^2)* M_w(i+1)*D2u(i,x);
end

