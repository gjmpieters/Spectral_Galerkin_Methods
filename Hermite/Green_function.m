function [y] = Green_function(eta, xi, b)

y = zeros(length(eta), length(xi));

for i=1:length(eta)
    for j=1:length(xi)
        if(eta(i) < xi(j))
            y(i,j) = (-1/(2*b))*(exp(b*(eta(i)-xi(j))) + exp(-b*(eta(i)+xi(j))));
        else
            y(i,j) = (-1/(2*b))*(exp(b*(xi(j)-eta(i))) + exp(-b*(xi(j)+eta(i))));
        end
    end
end


