function VsOut=L_SH(Vs,Ps,Es)
% Swift-Hohenberg model - Local terms
% VsOut=L_SH(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: U(1) : dU/dt = lambda*u + a2*u^2 + a3*u^3  + d2* D^2(u) + d4* D^4(u)
% Parameters are: lambda,a2,a3,d2,d4. (-1,2,-1,-2,-1)

if(~isfield(Es,'fmod'))
   Es.fmod=0;
end;

% Initialization
U=Vs(:,1); 

if(Es.fmod==0)      % Model equations

    dU = Ps.lambda.*U + Ps.a2.*U.*U + Ps.a3.*U.*U.*U;
    VsOut = dU;
else                % Jacobian of equations
    
    UdU = Ps.lambda + 2*Ps.a2.*U + 3*Ps.a3.*U.*U;
    % written in a large sparse matrix format 
    VsOut = sparse(diag(UdU));
end;

end