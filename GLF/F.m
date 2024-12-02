function [F] = F(s, no, nw, M)

nom = M*s.^nw;
denom = ((1-s).^no + M*s.^nw);

F =  nom ./ denom;

end
