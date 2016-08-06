function VsOut=L_NLL(Vs,Ps,Es,varargin)
% (negative) Lejeune-Lefever model - Local terms
% VsOut=L_NLL(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1)
% Parameters are: mu,Lambda,L. (1.02,1.2,0.2)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Initialization

B=Vs(:,1); 

if(Es.Jcob==0)
% Model equations

dB = (1 - Ps.mu).*B + (Ps.Lambda - 1).*B.*B - B.*B.*B;

VsOut = -dB;

else
% Jacobian of equations
BdB = 1 - Ps.mu + 2*(Ps.Lambda - 1).*B - 3*B.*B;

% written in a large sparse matrix format 
VsOut = sparse(diag(-BdB));
end;

