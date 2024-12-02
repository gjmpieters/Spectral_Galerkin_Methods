function dydt = ODE(eta, y, no, nw, M, lambda, Nc)

S = y(1);
U = y(2);

vFprime = F_prime(S, no, nw, M);
vG = G(S, no, nw, M, lambda);
vGprime = G_prime(S, no, nw, M, lambda);

dydt = [U ; ((-eta + 2*vFprime - 4*Nc*vG - 4*Nc*eta*vGprime*U)*U)/(4*Nc*eta*vG)];

% S = y(1);
% U = y(2);
% 
% vFprime = F_prime(S, no, nw, M);
% vG = G(S, no, nw, M, lambda);
% vGprime = G_prime(S, no, nw, M, lambda);
% 
% dydt = [U ; ((-eta + 2*vFprime - 4*Nc*vG - 4*Nc*eta*vGprime*U)*U)/(4*Nc*eta*vG)];

% Bernhard's formulation
%dydt = [U ; ((-eta + vFprime - 2*Nc*vG - 2*Nc*eta*vGprime*U)*U)/(2*Nc*eta*vG)];
end

