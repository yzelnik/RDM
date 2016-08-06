function VsOut=L_NVSG(Vs,Ps,Es,varargin)
% (Negative version of) Very simplified Gilad model - Local terms
% VsOut=L_NVSG(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1) and W(2). 
% Parameters are: P,eta,kappa,mu,ni,lambda,gamma,rho,DB,DW (1,3.5,1,1,3.333,3.333,3.333,0,1,1000)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Initialization

B=Vs(:,1); 
W=Vs(:,2); 

if(Es.Jcob==0)
% Model equations

dB = Ps.ni.*B.*W.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa) - Ps.mu.*B; 
dW = Ps.P - Ps.lambda.*W.*(1-Ps.rho.*B) - Ps.gamma.*B.*W.*(1 + Ps.eta.*B).^2;

VsOut = - [dB,dW];   % Negative change is here

else
% Jacobian of equations
BdB = (-Ps.kappa.* Ps.mu + (1 + B.* Ps.eta).* (Ps.kappa + B.* (-2 - 4* B.* Ps.eta + 3* Ps.eta.* Ps.kappa)).* Ps.ni.* W)./Ps.kappa;
BdW = Ps.ni.*B.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa);
WdB = - Ps.gamma.*W.*(1 + 4* B.* Ps.eta + 3* B.^2.* Ps.eta.^2) + Ps.rho.*Ps.lambda.*W;
WdW = - Ps.lambda.*(1-Ps.rho.*B) - Ps.gamma.*B.*(1 + Ps.eta.*B).^2 ;

% written in a large sparse matrix format 
VsOut = - sparse([diag(BdB) diag(BdW) ; diag(WdB) diag(WdW)]);    % Negative change is here
end;

