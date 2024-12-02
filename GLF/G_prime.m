function [Gprime] = G_prime(s, no, nw, M, lambda)

Gprime =    -F_prime(s, no, nw, M).*((1-s).^no).*Leverett_prime(s, lambda) +...
          no*F(s, no, nw, M).*((1-s).^(no-1)).*Leverett_prime(s, lambda) +...
            -F(s, no, nw, M).*((1-s).^no).*Leverett_primeprime(s, lambda);
 
%Gprime = 0;

end
