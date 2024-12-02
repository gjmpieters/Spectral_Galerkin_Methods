function [F, Jac] = Solve(Phi, quadrature_nodes, quadrature_weights, D, N, s, no, nw, M, lambda, Nc)

[mI, mE, mF, mG, mH, mJ_I, mJ_E, mJ_F, mJ_G, mJ_H] = Build_matrices_BL(Phi, quadrature_nodes, quadrature_weights, D, N, s, no, nw, M, lambda);

F = [sum(s); -mE + 2*mF + 4*Nc*mG] - eye(N+1,1);
%F = [sum(s); -0.5*mE - mH] - eye(N+1,1);



if nargout > 1 % gradient required
Jac = [ones(1,N); -mJ_E + 2*mJ_F + 4*Nc*mJ_G];
%Jac = [ones(1,N); -0.5*mJ_E - mJ_H];

end

