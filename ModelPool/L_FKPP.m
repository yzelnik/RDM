function VsOut=L_FKPP(Vs,Ps,Es,varargin)
% Fisher-Kolmogorov-Petrovsky-Piscounov equation - Local terms
% VsOut=L_FKPP(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: U(1) : dU/dt = r*u* (1-u) + d* D^2(u))
% Parameters are: r,d. (0.5,1)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
U=Vs(:,1); 

if(Es.JacMode==0)      % Model equations

    dU = Ps.r .* U .* (1-U);
    VsOut = dU;
else                % Jacobian of equations
    
    UdU = Ps.r .* (1-2*U);
    % written in a large sparse matrix format 
    VsOut = spdiags(UdU,0,Ps.Nx*Ps.Ny,Ps.Nx*Ps.Ny);
end;

end