function [u] = BuildSolution_BL(z, coefs, beta)

N = length(coefs);

u = Basis_BL(N, z, beta)*coefs;
