function VsOut=L_NISHI(Vs,Ps,Es,varargin)
% Nishiura model for gas discharge:
% du/dt = k2*u-u^3-k3*u-k4*v+k1+Du*u_xx
% dv/dt = (u-v+Dv*v_xx)/tau
% dw/dt = u-v+Dw*w_xx
%
% Taken from the paper "HETEROGENEITY-INDUCED SPOT DYNAMICS FOR A THREE-COMPONENT REACTION-DIFFUSION SYSTEM"
% Parameters used in paper: DU,DV,DW = [9*1e-5,1e-3,1e-2] , k1=k2=2, k3=1, k4=8.5 tau=40
% Ps_NISHI=struct('LocFunc',@L_NISHI,'SpaFunc',@S_LD,'k1',2,'k2',2,'k3',1,'k4',8.5,'tau',40,'VarNum',3,'Ds',[9*1e-5,1e-3,1e-2],'Lx',200,'Ly',1,'Nx',100,'Ny',1);

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Initialization

U=Vs(:,1); 
V=Vs(:,2); 
W=Vs(:,3); 

if(Es.JacMode==0)      % Model equations
    
    dU = Ps.k2.*U - U.^3 -Ps.k3.*V - Ps.k4.*W + Ps.k1;
    dV = (U-V)./Ps.tau;
    dW = U-W;
    
    VsOut = [dU,dV,dW];
else
% Jacobian of equations

% written in a large sparse matrix format 
VsOut = sparse(Ps.Nx*Ps.VarNum,Ps.Nx*Ps.VarNum);
warning ('Jacobian Not implemented.')
end;


