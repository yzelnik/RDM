function VsOut=L_EDGE(Vs,Ps,Es,varargin)
% Vegetation on the edge between UV and BS (based on the SG model) - Local terms
% That is, very simplified Gilad model with strong shading (dimensional), and extra B
% VsOut=L_EDGE(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1), W(2), H(3) and E(4)
% Ps = struct('LocFunc',@L_EDGE,'SpaFunc',@S_RD,'IntegFunc',@I_FDSIMP,'P',43,'eta',4.5,'kappa',0.5,'mu',0.25,'mu2',0.25,'nu',1.75,'lambda',0.005,'lambda2',0.005,'gamma',1.5,'gamma2',1.5,'rho',0.5,'Ds',[0.1 5 0.1],'Vnum',3,'Lx',160,'Ly',10,'Nx',160*5,'Ny',1,'BC',1);
% Parameters are: P,q,nu,alpha,eta,gamma,rho,f,DW,DH. (1,0.05 3.333 33.333 3.5 16.667  0.95 0.1,100,10000)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if(~isfield(Es,'fmod'))
   Es.fmod=0;
end;

% Initialization
B=Vs(:,1); 
W=Vs(:,2); 
H=Vs(:,3); 
E=Vs(:,4); 

if(Es.fmod==0)      % Model equations
    
    dB = Ps.lambda.*W.*B.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa) - Ps.mu.*B;
    dW = Ps.alpha.*(B*Ps.Y1 + E*Ps.Y2 + Ps.q.*Ps.f)./(B*Ps.Y1 + E*Ps.Y2 + Ps.q).*H - Ps.nu.*W./(1 + Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*W.*(1 + Ps.eta.*B).^2 - Ps.gamma2.*E.*W;
    dH = Ps.P - Ps.alpha.*(B*Ps.Y1 + E*Ps.Y2 + Ps.q.*Ps.f)./(B*Ps.Y1 + E*Ps.Y2 + Ps.q).*H;
    dE = Ps.lambda2.*E.*W.*(1 + Ps.eta2.*E).^2.*(1 - E./Ps.kappa2) - Ps.mu2.*E;
    
    VsOut = [dB,dW,dH,dE];
else               % Jacobian of equations
    BdB = (-Ps.kappa.* Ps.mu + (1 + B.* Ps.eta).* (Ps.kappa + B.* (-2 - 4* B.* Ps.eta + 3* Ps.eta.* Ps.kappa)).* Ps.lambda.* W)./Ps.kappa;
    BdW = Ps.lambda.*B.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa);
   % WdB = - Ps.gamma.*W.*(1 + 4* B.* Ps.eta + 3* B.^2.* Ps.eta.^2) + Ps.rho.*Ps.nu.*W./Ps.kappa;
   % WdW = - Ps.nu.*(1-Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*(1 + Ps.eta.*B).^2 ;
    % need to fix up rho part
    % written in a large sparse matrix format 
    VsOut = ArrangeJacobian([BdB BdW;WdB WdW],Ps,Es);
end;





end