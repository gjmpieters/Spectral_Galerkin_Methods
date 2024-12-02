function [Fprimeprime] = F_primeprime(s, no, nw, M)

n = ((1-s).^no + M*s.^nw).^2;
t = -M*(1-s).^(no-1).*s.^(nw-1).*((nw-no)*s - nw);

Fprimeprime = (n.*( M*(no-1) *(1-s).^(no-2).*s.^(nw-1).*((nw-no)*s - nw) +...
                   -M*(nw-1) *(1-s).^(no-1).*s.^(nw-2).*((nw-no)*s - nw) +...
                   -M*(nw-no)*(1-s).^(no-1).*s.^(nw-1) ) -...
               t.*( 2*(M*s.^nw + (1-s).^no).*(M*nw*s.^(nw-1) - no*(1-s).^(no-1)) ) ) ./ n.^2;

end
