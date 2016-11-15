function VsOut=L_GR(Vs,Ps,Es,varargin)
% ** WORK IN PROGRESS ** - Omer is trying to modify this one to Green Roofs model
% - local terms
% VsOut=L_OZFC(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B1(1),B2(2),S1(3). 
% Parameters are: P,q,mu,kappa,lambda,alpha,eta,gamma,nuW,rhoW,nuH,rhoH,f,DB,DW,DH. 
% Vals (non-dim): 2.5,3.6,1,1,1,13.333,1,1.554,0.5,0.3,1.5,0.8,0.01,1,25,2000
% Update online if necessary

% Ps=struct('LocFunc',@L_GR,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'P',80,'eta',6,'kappa',1,'mu',2,'nu',4,'lambda',0.06,'gamma',4.,'rho',0.2,'Ds',[0.02 20],'VarNum',2,'Lx',50,'Ly',1,'Nx',500,'Ny',1);
% Es=struct('TsSize',0.1,'TimeDst',200,'JacMode',0,'NonNeg',1,'LsaThresh',0,'OdeInit',1,'StSmall',0.001,'SsThresh',1e-5,'St1Color',circshift(hsv(6),-2)*0.66,'BfColor',[0 0 0 ; 1 0 0; 0 0.75 0; 0 0 1; 0 0.75 0.75; 1 0 1],'BfStyle',['-- ';'-  ';'x  ']);

if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
B1=Vs(:,1); 
B2=Vs(:,2); 
S1=Vs(:,3); 
if(Es.JacMode==0)      % Model equations

    dB1 = Ps.lamb1.*B1.*S1.*(1 + Ps.eta1.*B1).^2.*(1 - B1) - B1;
    dB2 = Ps.lamb2.*B2.*S1.*(1 + Ps.eta2.*B2).^2.*(1 - B2) - Ps.mu2.*B2;
    dS1 = Ps.p - Ps.nu.*S1./(1 + Ps.rho1*B1 + Ps.rho2*B2) - Ps.ks.*(S1.^4) - Ps.gam1.*B1.*S1.*(1 + Ps.eta1.*B1).^2 - Ps.gam2.*B2.*S1.*(1 + Ps.eta2.*B2).^2;

    VsOut = [dB1,dB2,dS1];
else                % Jacobian of equations
    
    B1dB1 = 2.*B1.*Ps.eta1.*Ps.lamb1.*S1.*(1 - B1).*(1 + Ps.eta1.*B1) - B1.*Ps.lamb1.*S1.*(1.0+B1.*Ps.eta1).^2 + Ps.lamb1.*S1.*(1 - B1).*(B1.*Ps.eta1 + 1).^2 - 1;
    B1dB2 = 0.*B1;
    B1dS1 = B1.*Ps.lamb1.*(1 - B1).*(B1.*Ps.eta1 + 1.0).^2;
    B2dB1 = 0.*B2;
    B2dB2 = 2*B2.*Ps.eta2.*Ps.lamb2.*S1.*(1 - B2).*(1 + Ps.eta2.*B2) - B2.*Ps.lamb2.*S1.*(B2.*Ps.eta2 + 1.0).^2 + Ps.lamb2.*S1.*(1 - B2).*(B2.*Ps.eta2 + 1.0).^2 - Ps.mu2;
    B2dS1 = B2.*Ps.lamb2.*(1 - B2).*(1 + Ps.eta1.*B2).^2;
    S1dB1 = -2*B1.*Ps.eta1.*Ps.gam1.*S1.*(1 + Ps.eta1.*B1) - Ps.gam1.*S1.*(1 + Ps.eta1.*B1).^2 + Ps.nu.*Ps.rho1.*S1./(B1.*Ps.rho1 + B2.*Ps.rho2 + 1.0).^2;
    S1dB2 = -2*B2.*Ps.eta2.*Ps.gam2.*S1.*(1 + Ps.eta2.*B2) - Ps.gam2.*S1.*(1 + Ps.eta2.*B2).^2 + Ps.nu.*Ps.rho2.*S1./(B1.*Ps.rho1 + B2.*Ps.rho2 + 1.0).^2;
    S1dS1 = -B1.*Ps.gam1.*(1 + Ps.eta1.*B1).^2 - B2.*Ps.gam2.*(1+ Ps.eta2.*B2).^2 - 4*Ps.ks.*(S1.^3) - Ps.nu./(B1.*Ps.rho1 + B2.*Ps.rho2 + 1.0);

    % written in a large sparse matrix format 
    VsOut = sparse([diag(B1dB1) diag(B1dB2) diag(B1dS1); diag(B2dB1) diag(B2dB2) diag(B2dS1); diag(S1dB1) diag(S1dB2) diag(S1dS1)]);
end;


