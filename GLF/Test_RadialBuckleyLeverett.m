%clear

% myCluster = parcluster('local');
% myCluster.NumWorkers = 8;  % 'Modified' property now TRUE
% saveProfile(myCluster);    % 'local' profile now updated,
                            % 'Modified' property now FALSE
% parpool(8);
                           
no = 4;
nw = 4;
M = 2;

lambda = 2 / (nw - 3);
%Nc = 15;
%Nc = 0.63;
Nc = 5;

% shooting method
%alpha = 1;
alpha = .5;
options = optimoptions('fsolve', 'Display', 'iter', 'MaxFunctionEvaluations', 15000, 'MaxIter',3000, 'StepTolerance', 1e-20, 'FunctionTolerance', 1e-20, 'OptimalityTolerance', 1e-20);
fun = @(alpha)SolveODE(alpha, no, nw, M, lambda, Nc);
[alphaopt, fval] = fsolve(fun, alpha, options)

options = odeset('Refine', 0, 'Stats', 'on', 'RelTol',1e-8,'AbsTol',1e-8);

fun = @(t, y)ODE(t, y, no, nw, M, lambda, Nc);

%eps = 1e-3;
%eta_start = alphaopt
%[t,y] = ode23(fun, [eta_start 0], [eps, -eps^(3/2)], options);
%[t,y] = ode23(fun, [eta_start eps], [eps, -0.00002], options);

%eta_start = (1-0.999)/-alphaopt;
%[t,y] = ode23(fun, [eta_start 10], [0.999, alphaopt], options);

eta_start = 0.125;
S_start = 1-alphaopt*eta_start^(1/no);
U_start = -(1/no)*alphaopt*eta_start^((1/no) - 1);
[t,y] = ode23s(fun, [eta_start 20], [S_start, U_start], options);

tt = [0:0.001:eta_start];

area = eta_start - (alphaopt/((1/no)+1))*eta_start^((1/no)+1) + sum(0.5*(y(2:end,1) + y(1:end-1,1)).*abs(t(2:end) - t(1:end-1)))

figure(1)
plot(t,y(:,1), tt, 1-alphaopt*tt.^(1/no))
axis square
hold on

figure(2)
%plot(y(:,1),y(:,2))
plot(t,y(:,2))
%axis([0 1 -5 0])
axis square
hold on


% N = 32;
% s = zeros(N,1);
% s(1) = 1;
% eps = [0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001]';
% 

% beta = 20;
% [quadrature_nodes, quadrature_weights] = Laguerre_Gauss_quadrature(N,0); 
% [Phi] = Basis_BL(N,quadrature_nodes/beta, beta);
% [D] = Derivative(Phi, quadrature_nodes/beta, N);        
 
% project solution onto basis

% for j=0 : N-1
%     %s(j+1) = sum( (1-erf(0.5*quadrature_nodes)).*Phi(:,j+1).*quadrature_weights);
%     s(j+1) = sum( interp1(t,y(:,1), quadrature_nodes/beta, 'pchip', 0).*Phi(:,j+1).*quadrature_weights);
% end

%figure(1)
%eta = 0:0.01:3;
%S = BuildSolution_BL(eta, s', beta);
%plot(eta, S);

%null(-mJ_E - mJ_H);

%options = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient',true, 'UseParallel', true, 'CheckGradients', false, 'StepTolerance', 1e-10, 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10)

%fun = @(s)Solve(Phi, quadrature_nodes, quadrature_weights, D, N, s, no, nw, M, lambda, Nc);
%[sopt, fval] = fsolve(fun, s, options)

% [mI, mE, mF, mG, mH, mJ_I, mJ_E, mJ_F, mJ_G, mJ_H] = Build_matrices_BL(Phi, N, s, no, nw, M, lambda);
% [mIm, mEm, mFm, mGm, mHm, mJ_Im, mJ_Em, mJ_Fm, mJ_Gm, mJ_Hm] = Build_matrices_BL(Phi, N,s - eps, no, nw, M, lambda);
% [mIp, mEp, mFp, mGp, mHp, mJ_Ip, mJ_Ep, mJ_Fp, mJ_Gp, mJ_Hp] = Build_matrices_BL(Phi, N,s + eps, no, nw, M, lambda);
% 
% mEp - mEm
% mJ_E*(2*eps)
% 
% 
% mFp - mFm
% mJ_F*(2*eps)
% 
% mGp - mGm
% mJ_G*(2*eps)

%eta = 0:0.1:70;
%S = BuildSolution_BL(eta, sopt);
%sum(0.5*(S(2:end,1) + S(1:end-1,1))'.*(eta(2:end) - eta(1:end-1)))

%figure(1)
%plot(quadrature_nodes, mB*sopt, '*');
%hold on;
%plot(eta, S);
%plot(eta,1-erf(eta/2))

