function VsOut=L_GS(Vs,Ps,Es)
% Gray-Scott model - Local terms
% VsOut=L_GS(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: V,U. Prameters are: f,k,D. Example values: (0.01,0.12,100)
% Equations: dV/dt = U*V^2 - (f+k)*V + V'',  dU/dt = -U*V^2 + f*(1-U) + D*U''

% Initilization
V=Vs(:,1); 
U=Vs(:,2); 

if(Es.JacMode==0)      % Model equations
    
    dV = U.*V.^2 - (Ps.k+Ps.f).*V;
    dU = -U.*V.^2 + Ps.f.*(1-U);
    
    VsOut = [dV,dU];
else    % Jacobian of equations
    
    VdV = 2*U.*V - (Ps.f+Ps.k);
    VdU = V.^2;
    UdV = - 2*U.*V;
    UdU = - V.^2 - Ps.f;
    
    % written in a large sparse matrix format
    VsOut = ArrangeJacobian([VdV VdU;UdV UdU],Ps,Es);
 end;


