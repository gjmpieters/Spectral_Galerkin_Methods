function k = permeability(z, beta)

% exp_up
%k = 1 - (1-beta)*exp(-beta * z);

% exp_down
%k = exp(-beta * z);

% step_down
%k = beta + (1-beta)*(1-heaviside(z-20));

% step_up
%k = beta + (1-beta)*heaviside(z-20);

% low
%k = beta;

% high
k = 1;

%%%% test %%%%k = 1 + 0.5*(beta-1) - 0.5*(beta-1)*tanh(0.5*(z-100));

 


