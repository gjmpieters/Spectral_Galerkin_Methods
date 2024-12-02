function F = SolveODE3(alpha, no, nw, M, lambda, Nc)

options = odeset('Refine', 0, 'Stats', 'on', 'RelTol',1e-3,'AbsTol',1e-3);

eta_start = (1-0.99)/-alpha(1);

fun = @(eta, y)ODE(eta, y, no, nw, M, lambda, Nc);
[t1,y1] = ode23(fun, [eta_start 2], [0.99, alpha(1)], options);


eta_start = alpha(2);

fun = @(eta, y)ODE(eta, y, no, nw, M, lambda, Nc);
[t2,y2] = ode23(fun, [eta_start 2], [0.001, -1000], options);

plot(t1,y1(:,1))
hold on
plot(t2,y2(:,1))

F = (2 - sum(0.5*(y1(2:end,1) + y1(1:end-1,1)).*abs(t1(2:end) - t1(1:end-1))) - sum(0.5*(y2(2:end,1) + y2(1:end-1,1)).*abs(t2(2:end) - t2(1:end-1))))^2 + (y1(end,1) - y2(end,1))^2 + (y1(end,2) - y2(end,2))^2;

%F = (y1(end,1) - y2(end,1))^2 + (y1(end,2) - y2(end,2))^2;

%F = (y1(end,1) - y2(end,1))^2;
