function [G] = G(s, no, nw, M, lambda)

G = -F(s, no, nw, M).*((1-s).^no).*Leverett_prime(s, lambda);

%G = 1e-2;

end
