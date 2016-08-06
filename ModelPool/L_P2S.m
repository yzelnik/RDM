function VsOut=L_P2S(Vs,Ps,Es,varargin)
% Paris's 2 species model - Local terms
% VsOut=L_P2S(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B1(1),B2(2),W(3). 
% Parameters are: P,ni,eta1,eta2,K1,K2,mu,L,gamma,DB1,DB2,DW (12,1,0.5,3,3,0.5,2,6,2.2,4.5,0.0000625,0.0000625,0.625)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Initialization

B1=Vs(:,1); 
B2=Vs(:,2); 
W =Vs(:,3); 

if(Es.fmod==0)      % Model equations

    dB1 = Ps.ni.*B1.*(1 + Ps.eta1.*B1).*(1 - B1./Ps.K1).*(1 - B2./(B1 + B2 + Ps.xhi)).*W - Ps.mu.*B1; 
    dB2 = Ps.ni.*B2.*(1 + Ps.eta2.*B2).*(1 - B2./Ps.K2).*(1 - B1./(B1 + B2 + Ps.xhi)).*W - Ps.mu.*B2; 
    dW  = Ps.P - Ps.L.*W - Ps.gamma.*W.*(B1.*(1 + Ps.eta1.*B1) + B2.*(1 + Ps.eta2.*B2));
    
    VsOut = [dB1,dB2,dW];
else                % Jacobian of equations

    B1dB1 = - Ps.mu + B1.*B2.*W.*(1 + B1.*Ps.eta1).*(1 - B1.*Ps.K1.^(-1)).*Ps.ni.*(B1 + B2 + Ps.xhi).^(-2) + B1.*W.*Ps.eta1.*(1 - B1.*Ps.K1.^(-1)).*Ps.ni.*(1 - B2.*(B1 + B2 + Ps.xhi).^(-1)) + W.*(1 + B1.*Ps.eta1).*(1 - B1.*Ps.K1.^(-1)).*Ps.ni.*(1 - B2.*(B1 + B2 + Ps.xhi).^(-1)) - B1.*W.*(1 + B1.*Ps.eta1).*Ps.K1.^(-1).*Ps.ni.*(1 - B2.*(B1 + B2 + Ps.xhi).^(-1));
    B1dB2 = B1.*W.*(1 + B1.*Ps.eta1).*(1 - B1.*Ps.K1.^(-1)).*Ps.ni.*(B2.*(B1 + B2 + Ps.xhi).^(-2) - (B1 + B2 + Ps.xhi).^(-1));
    B1dW  = B1.*(1 + B1.*Ps.eta1).*(1 - B1.*Ps.K1.^(-1)).*Ps.ni.*(1 - B2.*(B1 + B2 + Ps.xhi).^(-1));

    B2dB1 = B2.*W.*(1 + B2.*Ps.eta2).*(1 - B2.*Ps.K2.^(-1)).*Ps.ni.*(B1.*(B1 + B2 + Ps.xhi).^(-2) - (B1 + B2 + Ps.xhi).^(-1));
    B2dB2 = - Ps.mu + B1.*B2.*W.*(1 + B2.*Ps.eta2).*(1 - B2.*Ps.K2.^(-1)).*Ps.ni.*(B1 + B2 + Ps.xhi).^(-2) + B2.*W.*Ps.eta2.*(1 - B2.*Ps.K2.^(-1)).*Ps.ni.*(1 - B1.*(B1 + B2 + Ps.xhi).^(-1)) + W.*(1 + B2.*Ps.eta2).*(1 - B2.*Ps.K2.^(-1)).*Ps.ni.*(1 - B1.*(B1 + B2 + Ps.xhi).^(-1)) - B2.*W.*(1 + B2.*Ps.eta2).*Ps.K2.^(-1).*Ps.ni.*(1 - B1.*(B1 + B2 + Ps.xhi).^(-1));
    B2dW  = B2.*(1 + B2.*Ps.eta2).*(1 - B2.*Ps.K2.^(-1)).*Ps.ni.*(1 - B1.*(B1 + B2 + Ps.xhi).^(-1));

    WdB1  = - W.*(1 + 2.*B1.*Ps.eta1).*Ps.gamma;
    WdB2  = - W.*(1 + 2.*B2.*Ps.eta2).*Ps.gamma;
    WdW   = - (B1.*(1 + B1.*Ps.eta1) + B2.*(1 + B2.*Ps.eta2)).*Ps.gamma - Ps.L;

    % written in a large sparse matrix format 
    VsOut = sparse([diag(B1dB1) diag(B1dB2) diag(B1dW); diag(B2dB1) diag(B2dB2) diag(B2dW); diag(WdB1) diag(WdB2) diag(WdW)]);
end;






















