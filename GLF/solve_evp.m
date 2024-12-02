function [E,M] = solve_sym_evp(A,output,txt)

if output == 1
  txt = [' ',txt];
  T = strcat('Solving the ',txt,' eigenvalue problem...');
  disp(T)
end

N = length(A);

[Meig, Eeig] = eig(A,'balance');

%OPT.disp=0;
%[Meig Eeig] = eigs(A,N,'largestreal',OPT);

%OPT=struct('Disp',0,'LSolver','bicgstab');
%[Meig Eeig] = jdqr(inv(B)*A,1,'LR',OPT);

E = Eeig;
M = Meig;