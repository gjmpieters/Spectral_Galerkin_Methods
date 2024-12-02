function gs = groundstate(z,t,chi)

gs = (1/2)*erfc((z + chi*t)./(2*sqrt(t))) + (1/2)*exp(-chi*z).*erfc((z - chi * t)./(2*sqrt(t)));

gs(isnan(gs)) = 0;

% equilibrium
%gs = exp(-z);
