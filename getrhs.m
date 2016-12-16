function rhs=getrhs(Vs,Ps,Es,varargin)
% Get the Right-Hand-Side of the equation defined by Ps.LocFunc & Ps.SpaFunc
% rhs=getrhs(Vs,Ps,Es)

if(~mod(nargin,2)) error('No default extra-input exists for getrhs.'); end;

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Initilize state if necessary
[Vs,Ps,Es]=InitilizeState(Vs,Ps,Es);
% Setup the spatial matrix and auxiliary flags
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);		

Ps.SpaMat=zeros(Ps.Nx*Ps.VarNum);
% Calculate the right-hand-side
rhs=RightHandSide(Vs,Ps,Es);

end