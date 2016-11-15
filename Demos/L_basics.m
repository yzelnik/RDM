function Out=L_basics(Vs,Ps,Es)
% we return the rhs (right-hand-side) of the equation: u_t = r*u*(1-u)
Out = Ps.r * Vs .* (1-Vs);

end