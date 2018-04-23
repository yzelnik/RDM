function VsOut=L_EDGE(Vs,Ps,Es)
% Vegetation on the edge between UV and BS (based on the SG model) - Local terms
% That is, very simplified Gilad model with strong shading (dimensional), and extra B
% VsOut=L_EDGE(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1), W(2), H(3) and E(4)
% Ps = struct('LocFunc',@L_EDGE,'SpaFunc',@S_RD,'IntegFunc',@I_FDSIMP,'P',43,'eta',4.5,'kappa',0.5,'mu',0.25,'mu2',0.25,'nu',1.75,'lambda',0.005,'lambda2',0.005,'gamma',1.5,'gamma2',1.5,'rho',0.5,'Ds',[0.1 5 0.1],'VarNum',3,'Lx',160,'Ly',10,'Nx',160*5,'Ny',1,'BC',1);
% Parameters are: P,q,nu,alpha,eta,gamma,rho,f,DW,DH. (1,0.05 3.333 33.333 3.5 16.667  0.95 0.1,100,10000)

% Initialization
B=Vs(:,1); 
W=Vs(:,2); 
H=Vs(:,3); 
E=Vs(:,4); 

if(Es.JacMode==0)      % Model equations
    
    dB = Ps.lambda.*B.*W.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa) - Ps.mu.*B;
    dW = Ps.alpha.*(B*Ps.Y1 + E*Ps.Y2 + Ps.q.*Ps.f)./(B*Ps.Y1 + E*Ps.Y2 + Ps.q).*H - Ps.nu.*W./(1 + Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*W.*(1 + Ps.eta.*B).^2 - Ps.gamma2.*E.*W.*(1 + Ps.eta2.*E).^2;
    dH = Ps.P - Ps.alpha.*(B*Ps.Y1 + E*Ps.Y2 + Ps.q.*Ps.f)./(B*Ps.Y1 + E*Ps.Y2 + Ps.q).*H;
    dE = Ps.lambda2.*E.*W.*(1 + Ps.eta2.*E).^2.*(1 - E./Ps.kappa2) - Ps.mu2.*E;
    
    VsOut = [dB,dW,dH,dE];
else               % Jacobian of equations
    BdB = -Ps.mu + (1 + B.*Ps.eta).*(Ps.kappa+B.*(-2-4.*B.*Ps.eta + 3.*Ps.eta.*Ps.kappa)).*Ps.lambda.*W./Ps.kappa;
    BdW = B.*(1 + B.*Ps.eta).^2.*(1 - B./Ps.kappa).*Ps.lambda;
    
    WdB = -Ps.alpha.*(-1 + Ps.f).*H.*Ps.q.*Ps.Y1./(B*Ps.Y1 + E*Ps.Y2 + Ps.q).^2 + Ps.nu.*W.*Ps.rho./(Ps.kappa.*(1 + Ps.rho.*B./Ps.kappa).^2) - Ps.gamma.*W.*(1 + 4*Ps.eta.*B + 3*Ps.eta.^2.*B.^2);
    WdW = - Ps.nu./(1 + Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*(1 + Ps.eta.*B).^2 - Ps.gamma2.*E.*(1 + Ps.eta2.*E).^2; 
    WdH = Ps.alpha.*(B*Ps.Y1 + E*Ps.Y2 + Ps.q.*Ps.f)./(B*Ps.Y1 + E*Ps.Y2 + Ps.q);
    WdE = -Ps.alpha.*(-1 + Ps.f).*H.*Ps.q.*Ps.Y2./(B*Ps.Y1 + E*Ps.Y2 + Ps.q).^2 - Ps.gamma2.*W.*(1 + 4*Ps.eta2.*E + 3*Ps.eta2.^2.*E.^2);
    
    HdB = Ps.alpha.*(-1 + Ps.f).*H.*Ps.q.*Ps.Y1./(B*Ps.Y1 + E*Ps.Y2 + Ps.q).^2;
    HdH = -Ps.alpha.*(B*Ps.Y1 + E*Ps.Y2 + Ps.q.*Ps.f)./(B*Ps.Y1 + E*Ps.Y2 + Ps.q);
    HdE = Ps.alpha.*(-1 + Ps.f).*H.*Ps.q.*Ps.Y2./(B*Ps.Y1 + E*Ps.Y2 + Ps.q).^2;
    
    EdW = E.*(1 + E.*Ps.eta2).^2.*(1 - E./Ps.kappa2).*Ps.lambda2;
    EdE = -Ps.mu2 + (1 + E.*Ps.eta2).*(Ps.kappa2+E.*(-2-4.*E.*Ps.eta2 + 3.*Ps.eta2.*Ps.kappa2)).*Ps.lambda2.*W./Ps.kappa2;
    
    zrs = zeros(size(B));
  
    % written in a large sparse matrix format 
    VsOut = ArrangeJacobian([BdB BdW zrs zrs;WdB WdW WdH WdE;HdB zrs HdH HdE;zrs EdW zrs EdE],Ps,Es);
end;





end