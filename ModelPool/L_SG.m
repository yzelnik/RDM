function VsOut=L_SG(Vs,Ps,Es)
% Simplfied Gilad (non-dimensional) model - Local terms
% VsOut=L_SG(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1),W(2),H(3). 
% Parameters are: P,q,nu,alpha,eta,gamma,rho,f,DW,DH. (1,0.05 3.333 33.333 3.5 16.667  0.95 0.1,100,10000)

% Initialization
B=Vs(:,1); 
W=Vs(:,2); 
H=Vs(:,3); 

if(Es.JacMode==0)      % Model equations

    dB = Ps.nu.*W.*B.*(1 + Ps.eta.*B).^2.*(1 - B) - B;
    dW = Ps.alpha.*(B + Ps.q.*Ps.f)./(B + Ps.q).*H - Ps.nu*W.*(1 - Ps.rho*B) - Ps.gamma.*W.*(1 + Ps.eta.*B).^2.*B;
    dH = Ps.P - Ps.alpha.*(B + Ps.q.*Ps.f)./(B + Ps.q).*H;

    VsOut = [dB,dW,dH];
else                % Jacobian of equations
    
    BdB = -1-W.*(1+B.*Ps.eta).*(-1+2*B+B.*(-3+4*B).*Ps.eta).*Ps.nu;
    BdW = -(-1+B).*B.*(1+B.*Ps.eta).^2.*Ps.nu;
    BdH = 0.*B;
    WdB = -((-1+Ps.f).*H.*Ps.q.*Ps.alpha)./((B+Ps.q).^2)-W.*Ps.gamma.*(1+B.*Ps.eta).*(1+3*B.*Ps.eta)+Ps.rho.*W.*Ps.nu;
    WdW = -B.*Ps.gamma.*(1+B.*Ps.eta).^2+(-1+B.*Ps.rho).*Ps.nu;
    WdH = ((B+Ps.f.*Ps.q).*Ps.alpha)./(B+Ps.q);
    HdB = ((-1+Ps.f).*H.*Ps.q.*Ps.alpha)./((B+Ps.q).^2);
    HdW = 0.*H;
    HdH = -((B+Ps.f.*Ps.q).*Ps.alpha)./(B+Ps.q);

    % written in a large sparse matrix format 
    VsOut = ArrangeJacobian([BdB BdW BdH;WdB WdW WdH; HdB HdW HdH],Ps,Es);
end;


