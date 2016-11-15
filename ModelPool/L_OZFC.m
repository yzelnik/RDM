function VsOut=L_OZFC(Vs,Ps,Es,varargin)
% Simplfied Gilad (non-dimensional) for Aussie-Fairy-Circles model - Local terms
% VsOut=L_OZFC(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1),W(2),H(3). 
% Parameters are: P,q,mu,kappa,lambda,alpha,eta,gamma,nuW,rhoW,nuH,rhoH,f,DB,DW,DH. 
% Vals (non-dim): 2.5,3.6,1,1,1,13.333,1,1.554,0.5,0.3,1.5,0.8,0.01,1,25,2000
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

    dB = Ps.lambda.*W.*B.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa) - Ps.mu.*B;
    dW = Ps.alpha.*(B + Ps.q.*Ps.f)./(B + Ps.q).*H - Ps.nuW.*W./(1 + Ps.rhoW*B./Ps.kappa) - Ps.gamma.*W.*B.*(1 + Ps.eta.*B).^2;
    dH = Ps.P - Ps.alpha.*(B + Ps.q.*Ps.f)./(B + Ps.q).*H - Ps.nuH.*H./(1 + Ps.rhoH*B./Ps.kappa) ;

    VsOut = [dB,dW,dH];
else                % Jacobian of equations
    
    BdB = -Ps.mu + (1 + B.*Ps.eta).*(1 + 3*Ps.eta.*B + B./Ps.kappa.*(-2-4.*B.*Ps.eta) ).*Ps.lambda.*W;
    BdW = B.*(1 + B.*Ps.eta).^2.*(1 - B./Ps.kappa).*Ps.lambda;
    BdH = 0.*B;
    WdB = -Ps.alpha.*(Ps.f-1).*H.*Ps.q./(B+Ps.q).^2 + Ps.nuW.*Ps.rhoW./Ps.kappa.*W./(1 + Ps.rhoW*B./Ps.kappa).^2 - Ps.gamma.*W.*(1 + 4*Ps.eta.*B + 3*(Ps.eta*B).^2);
    WdW = - Ps.nuW./(1 + Ps.rhoW*B./Ps.kappa) - Ps.gamma.*B.*(1 + Ps.eta.*B).^2;
    WdH = Ps.alpha.*(B + Ps.f.*Ps.q)./(B + Ps.q);
    HdB = Ps.alpha.*(Ps.f-1).*H.*Ps.q./(B + Ps.q).^2 + Ps.nuH.*Ps.rhoH./Ps.kappa.*W./(1 + Ps.rhoH*B./Ps.kappa).^2;
    HdW = 0.*H;
    HdH = -Ps.alpha.*(B + Ps.f.*Ps.q)./(B + Ps.q) - Ps.nuH./(1 + Ps.rhoH*B./Ps.kappa);

    % written in a large sparse matrix format 
    VsOut = sparse([diag(BdB) diag(BdW) diag(BdH); diag(WdB) diag(WdW) diag(WdH); diag(HdB) diag(HdW) diag(HdH)]);
end;


