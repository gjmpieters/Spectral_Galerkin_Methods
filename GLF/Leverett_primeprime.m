function [Lprimeprime] = Leverett_primeprime(s, lambda)

Lprimeprime = (-1/lambda)*((-1/lambda)-1)*s.^((-1/lambda)-2);

end
