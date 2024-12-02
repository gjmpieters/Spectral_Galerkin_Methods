function k = permeability_prime(z, beta)

%k = beta * (1-beta) * exp(-beta * z);

%k = -beta * exp(-beta * z);

k = 1;

%k = -0.5*(beta-1)*0.5*sech(0.5*(z-100)).^2;
