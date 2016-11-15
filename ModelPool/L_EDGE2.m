function VsOut=L_EDGE2(Vs,Ps,Es,varargin)
% vegetation on the edge between UV and BS (based on the VSG model) - Local terms
% That is, very simplified Gilad model with strong shading (dimensional), and extra B
% VsOut=L_EDGE2(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1) and W(2) and E(3)
% Parameters for dimensional: P=100,eta=7,kappa=0.5,mu=0.25,nu=1.75,lambda=0.005,gamma=1.5,rho=0.5,DB=0.1,DW=5
% Parameters for non-dimensional: P=1.1,eta=3.5,kappa=(1),mu=(1),nu=7,lambda=(7),gamma=3,rho=0.5,DB=(1),DW=50

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
B=Vs(:,1); 
W=Vs(:,2); 
E=Vs(:,3);

if(Es.JacMode==0)      % Model equations

    dB = Ps.lambda.*B.*W.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa) - Ps.mu.*B; 
    dW = Ps.P                                      - Ps.nu.*W./(1 + Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*W.*(1 + Ps.eta.*B).^2 - Ps.gamma2.*E.*W;
    dE = Ps.lambda2.*E.*W - Ps.mu2.*E;
    
    VsOut = [dB,dW,dE];
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