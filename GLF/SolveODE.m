function F = SolveODE(alpha, no, nw, M, lambda, Nc)

options = odeset('Refine', 0, 'Stats', 'on', 'RelTol',1e-8,'AbsTol',1e-8);

eta_start = 0.125;
S_start = 1-alpha*eta_start^(1/no);
U_start = -(1/no)*alpha*eta_start^((1/no) - 1);

fun = @(eta, y)ODE(eta, y, no, nw, M, lambda, Nc);
[t,y] = ode23s(fun, [eta_start 20], [S_start, U_start], options);

area1 = eta_start - (alpha/((1/no)+1))*eta_start^((1/no)+1);
area2 = sum(0.5*(y(2:end,1) + y(1:end-1,1)).*abs(t(2:end) - t(1:end-1)));
F = y(end,1)^2 + (2 - (area1 + area2))^2;
