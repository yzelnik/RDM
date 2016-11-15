function jac=getjac(Vs,Ps,Es,varargin)
% Get the Right-Hand-Side of the equation defined by Ps.LocFunc & Ps.SpaFunc
% jac=getjac(Vs,Ps,Es)

% Default first extra input is for the Numerical-Jacobian flag
if(~mod(nargin,2)) varargin = ['Es.JacNum' varargin]; end;

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Initilize state if necessary
[Vs,Ps,Es]=InitilizeState(Vs,Ps,Es);
% Setup the spatial matrix and auxiliary flags
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);		

% Calculate the Jacobian
jac=CalculateJacobian(Vs,Ps,Es);

end