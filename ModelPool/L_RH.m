function VsOut=L_RH(Vs,Ps,Es)
% Rietkerk-HilleRisLambers (non-dimensional) model - Local terms
% VsOut=L_MK(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1),W(2),H(3). 
% Parameters are: P,mu,alpha,f,ni,gamma,DW,DH. (0.3,0.5,0.4,0.2,0.4,0.1,1,1000)

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
B=Vs(:,1); 
W=Vs(:,2); 
H=Vs(:,3); 

if(Es.JacMode==0)      % Model equations

    dB = B.*W./(1 + W) - Ps.mu.*B;
    dW = Ps.alpha.*H.*(B + Ps.f)./(B + 1) - Ps.ni.*W - Ps.gamma.*B.*W./(1 + W);
    dH = Ps.P - Ps.alpha.*H.*(B + Ps.f)./(B + 1);

    VsOut = [dB,dW,dH];
else                % Jacobian of equations
    BdB = W./(1 + W) - Ps.mu;
    BdW = B./((1 + W).^2);
    BdH = 0*B;
    WdB = Ps.alpha.*H.*(1 - Ps.f)./((B + 1).^2) - Ps.gamma.*W./(1 + W) ;
    WdW = - Ps.ni - Ps.gamma.*B./((1 + W).^2);
    WdH = Ps.alpha.*(B + Ps.f)./(B + 1);
    HdB = - Ps.alpha.*H.*(1 - Ps.f)./((B + 1).^2);
    HdW = 0*H;
    HdH = - Ps.alpha.*(B + Ps.f)./(B + 1);

    % written in a large sparse matrix format 
    VsOut = sparse([diag(BdB) diag(BdW) diag(BdH); diag(WdB) diag(WdW) diag(WdH); diag(HdB) diag(HdW) diag(HdH)]);
end;

end
