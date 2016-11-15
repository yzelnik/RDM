function jac=CalculateJacobian(Vs,Ps,Es)
% Calculate a Jacobian (either from an analytical form, or numerically
% jac=CalculateJacobian(Vs,Ps,Es,varargin)

if(Es.JacNum==1)
	jac = NumericJacobian(Vs,Ps,Es);
else
	Es.JacMode = 1;	% Request a jacobian
	jac = Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es);
end;

end
