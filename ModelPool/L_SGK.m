function VsOut=L_SGK(Vs,Ps,Es,varargin)
% Simplfied Gilad with K (non-dimensional) model - Local terms
% this version of the model used a different rescaling of the parameters
% such that the carrying capacity (K) remains
% VsOut=L_SGK(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1),W(2),H(3). 
% Parameters are: P,q,nu,alpha,eta,gamma,rho,f,DW,DH. (1,0.05 3.333 33.333 3.5 16.667  0.95 0.1,100,10000)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
B=Vs(:,1); 
W=Vs(:,2); 
H=Vs(:,3); 

if(Es.JacMode==0)      % Model equations

    dB = Ps.nu.*W.*B.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa) - B;
    dW = Ps.alpha.*(B + Ps.q.*Ps.f)./(B + Ps.q).*H - Ps.nu*W.*(1 - Ps.rho*B./Ps.kappa) - Ps.gamma.*W.*(1 + Ps.eta.*B).^2.*B;
    dH = Ps.P - Ps.alpha.*(B + Ps.q.*Ps.f)./(B + Ps.q).*H;

    VsOut = [dB,dW,dH];
else                % Jacobian of equations
    
    BdB = -1 + (1 + B.*Ps.eta).*(Ps.kappa+B.*(-2-4.*B.*Ps.eta + 3.*Ps.eta.*Ps.kappa)).*Ps.nu.*W./Ps.kappa;
    BdW = B.*(1 + B.*Ps.eta).^2.*(1 - B./Ps.kappa).*Ps.nu;
    BdH = 0.*B;
    WdB = -Ps.alpha.*(-1 + Ps.f).*H.*Ps.q./(B+Ps.q).^2 - Ps.nu.*((1 + B.*Ps.eta).*(1 + 3.*B.*Ps.eta).*Ps.kappa - Ps.rho).*W./Ps.kappa;
    WdW = Ps.nu.*(-1 - B.*(1 + B.*Ps.eta).^2 + B.*Ps.rho./Ps.kappa);
    WdH = Ps.alpha.*(B + Ps.f.*Ps.q)./(B + Ps.q);
    HdB = Ps.alpha.*(-1 + Ps.f).*H.*Ps.q./(B + Ps.q).^2;
    HdW = 0.*H;
    HdH = -Ps.alpha.*(B + Ps.f.*Ps.q)./(B + Ps.q);

% written in a large sparse matrix format 
VsOut = ArrangeJacobian([BdB BdW BdH;WdB WdW WdH; HdB HdW HdH],Ps,Es);

end;


