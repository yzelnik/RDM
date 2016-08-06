function VsOut=L_VSG(Vs,Ps,Es)
% Very simplified Gilad model (dimensional) - Local terms
% VsOut=L_VSG(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1) and W(2). 
% Parameters for dimensional: P=100,eta=7,kappa=0.5,mu=0.25,nu=1.75,lambda=0.005,gamma=1.5,rho=0.5,DB=0.1,DW=5
% Parameters for non-dimensional: P=1.1,eta=3.5,kappa=(1),mu=(1),nu=7,lambda=(7),gamma=3,rho=0.5,DB=(1),DW=50

if(~isfield(Es,'fmod'))
   Es.fmod=0;
end;

% Initialization
B=Vs(:,1); 
W=Vs(:,2); 

if(Es.fmod==0)      % Model equations

    dB = Ps.lambda.*B.*W.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa) - Ps.mu.*B; 
    dW = Ps.P - Ps.nu.*W.*(1-Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*W.*(1 + Ps.eta.*B).^2;

    VsOut = [dB,dW];
else               % Jacobian of equations
    BdB = (-Ps.kappa.* Ps.mu + (1 + B.* Ps.eta).* (Ps.kappa + B.* (-2 - 4* B.* Ps.eta + 3* Ps.eta.* Ps.kappa)).* Ps.lambda.* W)./Ps.kappa;
    BdW = Ps.lambda.*B.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa);
    WdB = - Ps.gamma.*W.*(1 + 4* B.* Ps.eta + 3* B.^2.* Ps.eta.^2) + Ps.rho.*Ps.nu.*W./Ps.kappa;
    WdW = - Ps.nu.*(1-Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*(1 + Ps.eta.*B).^2 ;

    % written in a large sparse matrix format 
    VsOut = ArrangeJacobian([BdB BdW;WdB WdW],Ps,Es);
end;





end