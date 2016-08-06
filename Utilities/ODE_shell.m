function dy = ODE_shell(~,y,Ps,Es)
% Make as shell over a Local function, so it can be run using ODE matlab functions
% dy = ODE_shell(t,y,Ps,Es)

dy = Ps.LocFunc(y(:)',Ps,Es)';

end
