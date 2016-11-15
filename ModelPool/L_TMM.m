function VsOut=L_TMM(Vs,Ps,Es,varargin)
% Thomas Mechanism model (1976) - Local terms
% VsOut=L_TMM(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: U,V
% Parameters are: {a,b,alpha,rho,K,D}={150,100,1.5,13,0.05,100}

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Initialization

U=Vs(:,1); 
V=Vs(:,2); 

if(Es.JacMode==0)
% Model equations

dU = Ps.a - U - Ps.rho.*U.*V./(1 + U + Ps.K.*U.^2);
dV = Ps.alpha .*(Ps.b - V) - Ps.rho.*U.*V./(1 + U + Ps.K.*U.^2);

VsOut = [dU,dV];

else
% Jacobian of equations

% written in a large sparse matrix format 
VsOut = sparse([zeros(size(Vs))]);
end;


