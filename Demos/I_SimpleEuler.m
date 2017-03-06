function Vs=I_SimpleEuler(Vs,Ps,Es)
% Integrator with Finite-Difference Euler scheme

% How many time-steps are we running in our "for" loop?
totsteps = ceil(Es.TimeDst/Es.TsSize);

% Go through each time step
for ii=1:totsteps
    loc_rhs = Ps.LocFunc(Vs,Ps,Es);
    spa_rhs = Ps.SpaFunc(Vs,Ps,Es)*Vs;
    VsNew = Vs + Es.TsSize*(loc_rhs + spa_rhs);
    Vs = VsNew;
end;

end
