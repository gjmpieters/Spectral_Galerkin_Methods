function F = SolveODE2(alpha, no, nw, M, lambda, Nc)

options = odeset('Refine', 0, 'Stats', 'on', 'RelTol',1e-3,'AbsTol',1e-3);

eps = 1e-3;
eta_start = alpha;


fun = @(eta, y)ODE(eta, y, no, nw, M, lambda, Nc);
[t,y] = ode23(fun, [eta_start 0], [eps, -eps^(3/2)], options);
%[t,y] = ode23(fun, [eta_start eps], [eps, -0.00002], options);


%plot(t,y(:,1))
%hold on

%F = (y(end,1) - 1)^2

area = 2 - sum(0.5*(y(2:end,1) + y(1:end-1,1)).*abs(t(2:end) - t(1:end-1)));
F = (area)^2

%F = 2 + 0.001*eta_start - sum(0.5*(y(2:end,1) + y(1:end-1,1)).*abs(t(2:end) - t(1:end-1)))
