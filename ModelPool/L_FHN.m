function VsOut=L_FHN(Vs,Ps,Es,varargin)
% FitzHugh-Nagumo - Local terms
% VsOut=L_FHN(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: U,V
% Parameters are: {a0,a1,epsilon,delta}={-0.1,2,0.05,4}
% Ps_FHN=struct('LocFunc',@L_FHN,'SpaFunc',@S_RD,'a0',-0.1,'a1',2,'epsilon',0.05,'VarNum',2,'Ds',[1 4],'Lx',100,'Ly',100,'Nx',512,'Ny',1);

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Initialization
U=Vs(:,1); 
V=Vs(:,2); 

if(Es.JacMode==0)      % Model equations

    dU = U - U.^3 - V;
    dV = Ps.epsilon .*(U - Ps.a1.*V - Ps.a0);

    VsOut = [dU,dV];
else                % Jacobian of equations
    UdU = 1 - 2*U.^2;
    UdV = -1;
    VdU = Ps.epsilon;
    VdV = -Ps.epsilon.*Ps.a1;
    
    syslen = Ps.Nx * Ps.Ny;
    % written in a large sparse matrix format
    VsOut = spdiags([UdU UdV VdU; VdV UdV VdU],[0 syslen -syslen],syslen*2,syslen*2);
    %VsOut = sparse([diag(UdU) diag(UdV) ; diag(VdU) diag(VdV)]);
end;

end
