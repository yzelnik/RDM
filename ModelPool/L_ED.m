function VsOut=L_ED(Vs,Ps,Es,varargin)
% Exponential decay - Local terms
% VsOut=L_ED(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Parameter is lambda for half-life

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;


if(Es.JacMode==0)      % Model equations
    VsOut = -Vs./Ps.lambda;
else               % Jacobian of equations
    % written in a large sparse matrix format 
    VsOut = spdiags(repmat(-1./Ps.lambda,Ps.Nx*Ps.Ny,1),0,Ps.Nx*Ps.Ny,Ps.Nx*Ps.Ny);
end;


end