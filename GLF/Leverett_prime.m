function [Lprime] = Leverett_prime(s, lambda)

Lprime = (-1/lambda)*s.^((-1/lambda)-1);

end
