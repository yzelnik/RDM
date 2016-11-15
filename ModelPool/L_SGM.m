function VsOut=L_SGM(Vs,Ps,Es,varargin)
% Simplfied Gilad with Mu model - Local terms
% VsOut=L_SG(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1),W(2),H(3). 
% Parameters are: P,q,ni,alpha,eta,gamma,rho,f,DW,DH. (1,0.05 3.333 33.333 3.5 16.667  0.95 0.1,100,10000)
% Also, mu, default is 1

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Initialization

B=Vs(:,1); 
W=Vs(:,2); 
H=Vs(:,3); 

if(Es.JacMode==0)
% Model equations

dB = Ps.ni.*W.*B.*(1 + Ps.eta.*B).^2.*(1 - B) - Ps.mu*B;
dW = Ps.alpha.*(B + Ps.q.*Ps.f)./(B + Ps.q).*H - Ps.ni*W.*(1 - Ps.rho*B) - Ps.gamma.*W.*(1 + Ps.eta.*B).^2.*B;
dH = Ps.P - Ps.alpha.*(B + Ps.q.*Ps.f)./(B + Ps.q).*H;

VsOut = [dB,dW,dH];

else
% Jacobian of equations
BdB = -Ps.mu-W.*(1+B.*Ps.eta).*(-1+2*B+B.*(-3+4*B).*Ps.eta).*Ps.ni;
BdW = -(-1+B).*B.*(1+B.*Ps.eta).^2.*Ps.ni;
BdH = 0.*B;
WdB = -((-1+Ps.f).*H.*Ps.q.*Ps.alpha)./((B+Ps.q).^2)-W.*Ps.gamma.*(1+B.*Ps.eta).*(1+3*B.*Ps.eta)+Ps.rho.*W.*Ps.ni;
WdW = -B.*Ps.gamma.*(1+B.*Ps.eta).^2+(-1+B.*Ps.rho).*Ps.ni;
WdH = ((B+Ps.f.*Ps.q).*Ps.alpha)./(B+Ps.q);
HdB = ((-1+Ps.f).*H.*Ps.q.*Ps.alpha)./((B+Ps.q).^2);
HdW = 0.*H;
HdH = -((B+Ps.f.*Ps.q).*Ps.alpha)./(B+Ps.q);

% written in a large sparse matrix format 
VsOut = sparse([diag(BdB) diag(BdW) diag(BdH); diag(WdB) diag(WdW) diag(WdH); diag(HdB) diag(HdW) diag(HdH)]);
end;


