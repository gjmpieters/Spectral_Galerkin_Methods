function gs = primprof(z,t,chi)

if chi == -1 
% downward flow
gs = (1/2)*exp(z).*erfc((z+t)/(2*sqrt(t))) - (1/2)*exp(z).*exp(-(z+t).^2/(4*t))/(sqrt(pi*t)) - (1/2)*exp(-(z-t).^2/(4*t))/sqrt(pi*t);
gs = -gs;

elseif chi == 1 
% upward flow
gs = -(1/2)*exp(-z).*erfc((z-t)/(2*sqrt(t))) - (1/2)*exp(-z).*exp(-(z-t).^2/(4*t))/(sqrt(pi*t)) - (1/2)*exp(-(z+t).^2/(4*t))/sqrt(pi*t);
gs = - gs;

elseif chi == 0 
% no flow
gs = -exp(-z.^2/(4*t))/(sqrt(pi*t));
gs = -gs;
end

gs(isnan(gs)) = 0;

% equilibrium
%gs = exp(-z);


