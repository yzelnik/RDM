function VsOut=L_DCn(Vs,Ps,Es,varargin)
% Dynamic Crust (dimensional) model - Local terms
% VsOut=L_SG(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1),C(2),W(3),H(4). 
% Parameters are: {P,LambdaB,KB,MB,PhiB,LambdaC,KC,MC,PhiC,A,f,Q,N,GammaB,GammaC,DB,DC,DW,DH}={170,0.032,10,1.2,0.1,0.032,0.1,0.6,10,40,0.1,0.05,4,20,0.032,6.25*10^(-4),6.25*10^(-3),6.25*10^(-2),5*10^(0)}

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Initialization

B=Vs(:,1); 
C=Vs(:,2); 
W=Vs(:,3); 
H=Vs(:,4); 

if(Es.JacMode==0)
% Model equations
II = Ps.alpha.*(Ps.qc+C.*Ps.fc)./(Ps.qc+C).*((B+Ps.qb*Ps.fb)./(B+Ps.qb));

dB = Ps.nu.*W.*(1+Ps.eta.*B).*(1+Ps.eta.*B).*B.*(1 - B) - B - (Ps.phib.*Ps.k.*B.*(C.^Ps.nn))./((B+Ps.bZero).^Ps.mm);
dC = Ps.nu.*(Ps.lambdacw.*W+Ps.lambdach.*H).*C.*(1 - C) - Ps.mu.*C - Ps.phic.*B.*C;
dW = II.*H - Ps.nu.*(1-Ps.r.*B).*W - (Ps.gammab.*B.*(1+Ps.eta.*B).*(1+Ps.eta.*B) + Ps.gammacw.*C).*W;
dH = Ps.p - II.*H - Ps.gammach.*C.*H;

VsOut = [dB,dC,dW,dH];

else
% Jacobian of equations
BdB = -1-W.*(1+B.*Ps.eta).*(-1+2.*B+B.*(-3+4.*B).*Ps.eta).*Ps.nu-C.^(Ps.nn).*(B+Ps.bZero).^(-1-Ps.mm).*Ps.k.*(B+Ps.bZero-B.*Ps.mm).*Ps.phib;
BdC = -B.*C.^(-1+Ps.nn).*(B+Ps.bZero).^(-Ps.mm).*Ps.k.*Ps.nn.*Ps.phib;
BdW = -(-1+B).*B.*(1+B.*Ps.eta).^2.*Ps.nu;
BdH = 0.*B;
CdB = -C.*Ps.phic;
CdC = -Ps.mu-(-1+2.*C).*(H.*Ps.lambdach+W.*Ps.lambdacw).*Ps.nu-B.*Ps.phic;
CdW = -(-1+C).*C.*Ps.lambdacw.*Ps.nu;
CdH = -(-1+C).*C.*Ps.lambdach.*Ps.nu;
WdB = -W.*(1+B.*Ps.eta).*(1+3.*B.*Ps.eta).*Ps.gammab-(H.*Ps.alpha.*(-1+Ps.fb).*Ps.qb.*(C.*Ps.fc+Ps.qc))./((B+Ps.qb).^2.*(C+Ps.qc))+W.*Ps.nu.*Ps.r;
WdC = (H.*Ps.alpha.*(-1+Ps.fc).*(B+Ps.fb.*Ps.qb).*Ps.qc-W.*Ps.gammacw.*(B+Ps.qb).*(C+Ps.qc).^2)./((B+Ps.qb).*(C+Ps.qc).^2);
WdW = -B.*(1+B.*Ps.eta).^2.*Ps.gammab-C.*Ps.gammacw+Ps.nu.*(-1+B.*Ps.r);
WdH = (Ps.alpha.*(B+Ps.fb.*Ps.qb).*(C.*Ps.fc+Ps.qc))./((B+Ps.qb).*(C+Ps.qc));
HdB = (H.*Ps.alpha.*(-1+Ps.fb).*Ps.qb.*(C.*Ps.fc+Ps.qc))./((B+Ps.qb).^2.*(C+Ps.qc));
HdC = -((H.*(Ps.alpha.*(-1+Ps.fc).*(B+Ps.fb.*Ps.qb).*Ps.qc+Ps.gammach.*(B+Ps.qb).*(C+Ps.qc).^2))./((B+Ps.qb).*(C+Ps.qc).^2));
HdW = 0.*H;
HdH = -C.*Ps.gammach-(Ps.alpha.*(B+Ps.fb.*Ps.qb).*(C.*Ps.fc+Ps.qc))./((B+Ps.qb).*(C+Ps.qc));

%disp('not implemented yet');

% written in a large sparse matrix format 
VsOut = sparse([diag(BdB) diag(BdC) diag(BdW) diag(BdH); diag(CdB) diag(CdC) diag(CdW) diag(CdH); diag(WdB) diag(WdC) diag(WdW) diag(WdH); diag(HdB) diag(HdC) diag(HdW) diag(HdH)]);
end;


