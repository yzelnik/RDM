function VsOut=L_DC(Vs,Ps,Es,varargin)
% Dynamic Crust (dimensional) model - Local terms
% VsOut=L_SG(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1),C(2),W(3),H(4). 
% Parameters are: {P,LambdaB,KB,MB,PhiB,LambdaC,KC,MC,PhiC,A,f,Q,N,GammaB,GammaC,DB,DC,DW,DH}={170,0.032,10,1.2,0.1,0.032,0.1,0.6,10,40,0.1,0.05,4,20,0.032,6.25*10^(-4),6.25*10^(-3),6.25*10^(-2),5*10^(0)}

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Initialization

B=Vs(:,1); 
C=Vs(:,2); 
W=Vs(:,3); 
H=Vs(:,4); 

if(Es.fmod==0)
% Model equations
%II = Ps.A.*(Ps.QC+C.*Ps.fC)./(Ps.QC+C);
II = Ps.A.*(Ps.QC+C.*Ps.fC)./(Ps.QC+C).*((B+Ps.QB*Ps.fB)./(B+Ps.QB));

dB = Ps.LambdaB.*W.*B.*(1 - B./Ps.KB).*(1+Ps.EE.*B).*(1+Ps.EE.*B) - Ps.MB.*B - (Ps.PhiB.*B.*C)./(Ps.beta.*B+Ps.BZero);
dC = (Ps.LambdaCW.*W+Ps.LambdaCH.*H).*C.*(1 - C./Ps.KC) - Ps.MC.*C - Ps.PhiC.*B.*C;
dW = II.*H - Ps.N.*(1-Ps.R.*B./Ps.KB).*W - (Ps.GammaB.*B.*(1+Ps.EE.*B).*(1+Ps.EE.*B) + Ps.GammaCW.*C).*W;
dH = Ps.P - II.*H - Ps.GammaCH.*C.*H;

VsOut = [dB,dC,dW,dH];

else
% Jacobian of equations
BdB = -Ps.MB - (Ps.BZero.*C.*Ps.PhiB)./(B+Ps.BZero).^2;
BdC = -B.*Ps.PhiB;
BdW = B.*(1-B./Ps.KB).*Ps.LambdaB;
BdH = 0.*B;
CdB = -C.*Ps.PhiC;
CdC= -Ps.MC -B.*Ps.PhiC + (1-C./Ps.KC).*Ps.LambdaC.*W - (C.*Ps.LambdaC.*W)./Ps.KC;
CdW = C.*(1-C./Ps.KC).*Ps.LambdaC;
CdH = 0.*C;
WdB = -Ps.GammaC.*W;
WdC = -(Ps.A.*Ps.f.*H.*(C+Ps.Q))./(C+Ps.f.*Ps.Q).^2 + (Ps.A.*Ps.f.*H)./(C+Ps.f.*Ps.Q) - Ps.GammaC.*W;
WdW = -B.*Ps.GammaB - C.*Ps.GammaC - Ps.N;
WdH = (Ps.A.*Ps.f.*(C+Ps.Q))./(C+Ps.f.*Ps.Q);
HdB = 0.*H;
HdC = (Ps.A.*Ps.f.*H.*(C+Ps.Q))./(C+Ps.f.*Ps.Q).^2 - (Ps.A.*Ps.f.*H)./(C+Ps.f.*Ps.Q);
HdW = 0.*H;
HdH = -(Ps.A.*Ps.f.*(C+Ps.Q))./(C+Ps.f.*Ps.Q);

warning('Jacobian not implemented yet');

% written in a large sparse matrix format 
VsOut = sparse([diag(BdB) diag(BdC) diag(BdW) diag(BdH); diag(CdB) diag(CdC) diag(CdW) diag(CdH); diag(WdB) diag(WdC) diag(WdW) diag(WdH); diag(HdB) diag(HdC) diag(HdW) diag(HdH)]);
end;


