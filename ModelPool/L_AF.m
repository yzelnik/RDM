function VsOut=L_AF(Vs,Ps,Es,varargin)
% ** WORK IN PROGRESS ** - Omer is trying to modify this one to AgroForestry model
% - local terms
% VsOut=L_OZFC(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B1(1),B2(2),S1(3). 
% Parameters are: P,q,mu,kappa,lambda,alpha,eta,gamma,nuW,rhoW,nuH,rhoH,f,DB,DW,DH. 
% Vals (non-dim): 2.5,3.6,1,1,1,13.333,1,1.554,0.5,0.3,1.5,0.8,0.01,1,25,2000
% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'fmod'))
   Es.fmod=0;
end;

% Initialization
B1=Vs(:,1); 
B2=Vs(:,2); 
S1=Vs(:,3);
S2=Vs(:,4);
if(Es.fmod==0)      % Model equations

    dB1 = Ps.lamb1.*(1-(B2.^Ps.l_comp)./((Ps.k.*B1).^Ps.l_comp+B2.^Ps.l_comp+Ps.Br.^Ps.l_comp)).*B1.*S1.*(1 + Ps.eta1.*B1).^2.*(1 - B1) - B1;
    dB2 = Ps.lamb2.*(1-((Ps.k.*B1).^Ps.l_comp)./((Ps.k.*B1).^Ps.l_comp+B2.^Ps.l_comp+Ps.Br.^Ps.l_comp)).*B2.*S2.*(1 + Ps.eta2.*B2).^2.*(1 - B2) - Ps.mu2.*B2;
    dS1 = Ps.p - Ps.nu.*S1./(1 + Ps.rho1*B1 + Ps.rho2*B2) - Ps.ks.*(S1.^4) - Ps.gam1.*B1.*S1.*(1 + Ps.eta1.*B1).^2;
    dS2 = (Ps.ks./Ps.zeta)*((S1.^4) - (S2.^4)) - Ps.gam2.*B2.*S2.*(1 + Ps.eta2.*B2).^2;
    
    VsOut = [dB1,dB2,dS1,dS2];
else                % Jacobian of equations
    %B1dB1 =  2*b1*eta1*lam1*s1*(-b1 + 1.0)*(b1*eta1 +
    %1.0)*(-b2**l_comp/(b1**l_comp*k + b2**l_comp + br**l_comp) + 1.0) - b1*lam1*s1*(b1*eta1 + 1.0)**2*(-b2**l_comp/(b1**l_comp*k + b2**l_comp + br**l_comp) + 1.0) + b1**l_comp*b2**l_comp*k*l_comp*lam1*s1*(-b1 + 1.0)*(b1*eta1 + 1.0)**2/(b1**l_comp*k + b2**l_comp + br**l_comp)**2 + lam1*s1*(-b1 + 1.0)*(b1*eta1 + 1.0)**2*(-b2**l_comp/(b1**l_comp*k + b2**l_comp + br**l_comp) + 1.0) - 1;
    %B1dB2 =  b1*lam1*s1*(-b1 + 1.0)*(b1*eta1 +
    %1.0)**2*(b2**(2*l_comp)*l_comp/(b2*(b1**l_comp*k + b2**l_comp + br**l_comp)**2) - b2**l_comp*l_comp/(b2*(b1**l_comp*k + b2**l_comp + br**l_comp)));
    %B1dS1 =  b1*lam1*(-b1 + 1.0)*(b1*eta1 +
    %1.0)**2*(-b2**l_comp/(b1**l_comp*k + b2**l_comp + br**l_comp) + 1.0);
    %B1dS2 =  0;
    %B2dB1 =  b2*lam2*s2*(-b2 + 1.0)*(b2*eta2 +
    %1.0)**2*(b1**(2*l_comp)*k**2*l_comp/(b1*(b1**l_comp*k + b2**l_comp + br**l_comp)**2) - b1**l_comp*k*l_comp/(b1*(b1**l_comp*k + b2**l_comp + br**l_comp)));
    %B2dB2 =  b1**l_comp*b2**l_comp*k*l_comp*lam2*s2*(-b2 + 1.0)*(b2*eta2 +
    %1.0)**2/(b1**l_comp*k + b2**l_comp + br**l_comp)**2 + 2*b2*eta2*lam2*s2*(-b2 + 1.0)*(b2*eta2 + 1.0)*(-b1**l_comp*k/(b1**l_comp*k + b2**l_comp + br**l_comp) + 1.0) - b2*lam2*s2*(b2*eta2 + 1.0)**2*(-b1**l_comp*k/(b1**l_comp*k + b2**l_comp + br**l_comp) + 1.0) + lam2*s2*(-b2 + 1.0)*(b2*eta2 + 1.0)**2*(-b1**l_comp*k/(b1**l_comp*k + b2**l_comp + br**l_comp) + 1.0) - mu2;
    %B2dS1 =  0;
    %B2dS2 =  b2*lam2*(-b2 + 1.0)*(b2*eta2 +
    %1.0)**2*(-b1**l_comp*k/(b1**l_comp*k + b2**l_comp + br**l_comp) + 1.0);
    B1dB1 = 2.*B1.*Ps.eta1.*Ps.lamb1.*S1.*(1 - B1).*(1 + Ps.eta1.*B1) - B1.*Ps.lamb1.*S1.*(1 + Ps.eta1.*B1).^2 + Ps.lamb1.*S1.*(1 - B1).*(1 + Ps.eta1.*B1).^2 - 1;
    B1dB2 = 0.*B1;
    B1dS1 = B1.*Ps.lamb1.*(1 - B1).*(1 + Ps.eta1.*B1).^2;
    B1dS2 = 0.*B1;
    B2dB1 = 0.*B2;
    B2dB2 = 2.*B2.*Ps.eta2.*Ps.lamb2.*S2.*(1 - B2).*(1 + Ps.eta2.*B2) - B2.*Ps.lamb2.*S2.*(1+ Ps.eta2.*B2).^2 + Ps.lamb2.*S2.*(1 - B2).*(1 + Ps.eta2.*B2).^2 - Ps.mu2;
    B2dS1 = 0.*B2;
    B2dS2 = B2.*Ps.lamb2.*(1 - B2).*(1 + Ps.eta2.*B2).^2;
    S1dB1 = -2.*B1*Ps.eta1.*Ps.gam1.*S1.*(1 + Ps.eta1.*B1) - Ps.gam1.*S1.*(1 + Ps.eta1.*B1).^2 + Ps.nu.*Ps.rho1.*S1;
    S1dB2 = Ps.nu.*Ps.rho2.*S1;
    S1dS1 = -B1.*Ps.gam1.*(1 + Ps.eta1.*B1).^2 - 4.*Ps.ks.*(S1.^3) - Ps.nu.*(1 -Ps.rho1.*B1 - Ps.rho2.*B2);
    S1dS2 = 0.*S1;
    S2dB1 = 0.*S2;
    S2dB2 = -2.*B2.*Ps.eta2.*Ps.gam2.*S2.*(1 + Ps.eta2.*B2) - Ps.gam2.*S2.*(1 + Ps.eta2.*B2).^2;
    S2dS1 = 4.*Ps.ks.*(S1.^3)./Ps.zeta;
    S2dS2 = -B2.*Ps.gam2.*(1 + Ps.eta2.*B2).^2 - 4.*Ps.ks.*(S2.^3)./zeta;

    % written in a large sparse matrix format 
    VsOut = sparse([diag(B1dB1) diag(B1dB2) diag(B1dS1) diag(B1dS2); diag(B2dB1) diag(B2dB2) diag(B2dS1) diag(B2dS2);
                    diag(S1dB1) diag(S1dB2) diag(S1dS1) diag(S1dS2); diag(S2dB1) diag(S2dB2) diag(S2dS1) diag(S2dS2)]);
end;


