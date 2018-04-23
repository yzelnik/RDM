function VsOut=L_CSV(Vs,Ps,Es)
% Coastal Salinity Vegetation model - Local terms
% VsOut=L_CSV(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), 
% calculate the local terms of the model
% Variables are: N1(1),N2(2),S(3). 
% Parameters are: p1,p2,a11,a12,a21,a22,b0,b1,m,k,e,g. 
% (0.1,0.1,0.03,0.02,0.06,0.03,0.3,1,3.14,1.2,1.5,5)

% Initialization

N1 = Vs(:,1); 
N2 = Vs(:,2); 
S  = Vs(:,3); 

if(Es.JacMode==0)
% Model equations

dN1 = N1.*(Ps.p1.*(Ps.m./(Ps.m + S)) - Ps.a11.*N1 - Ps.a12.*N2);
dN2 = N2.*(Ps.p2 - Ps.a21.*N1 - Ps.a22.*N2);
dS  = Ps.b0.*Ps.g + N2.*Ps.b1.*Ps.g./(Ps.k + N2) - Ps.e.*S;

VsOut = [dN1,dN2,dS];

else
% Jacobian of equations
warning('Jacobian not implemented');
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


