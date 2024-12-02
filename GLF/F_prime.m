function [Fprime] = F_prime(s, no, nw, M)

nom = M*s.^nw;
denom = (1-s).^no + M*s.^nw;

t = nw*M*s.^(nw-1) * denom - (-no*(1-s).^(no-1) + nw*M*s.^(nw-1)) .* nom;

Fprime = t ./ denom.^2;

end
