function VsOut=L_NDX(Vs,Ps,Es,varargin)
% Neutron Diffusion with Xenon
% VsOut=L_NDX(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: P(1),I(2),X(3). 
% Parameters are: P,q,ni,alpha,eta,gamma,rho,f,DW,DH. (1,0.05 3.333 33.333 3.5 16.667  0.95 0.1,100,10000)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
P=Vs(:,1); 
I=Vs(:,2); 
X=Vs(:,3); 

if(Es.JacMode==0)      % Model equations

    dP = Ps.v.*(Ps.nu.*Ps.Sigma_f.*P - Ps.Sigma_a.*P - Ps.sigma_x.*X.*P);
    dI = Ps.gamma_I.*Ps.Sigma_f.*P - Ps.lambda_I.*I;
    dX = Ps.gamma_X.*Ps.Sigma_f.*P + Ps.lambda_I.*I - Ps.lambda_X.*X - Ps.sigma_x.*X.*P;

    VsOut = [dP,dI,dX];
else                % Jacobian of equations
    
    BdB = -1-W.*(1+B.*Ps.eta).*(-1+2*B+B.*(-3+4*B).*Ps.eta).*Ps.ni;
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


