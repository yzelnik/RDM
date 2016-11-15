function VsOut=L_VSGE1(Vs,Ps,Es)
% Very simplified Gilad model, but with only (1+EB)^1
% VsOut=L_VSGE1(Vs,Ps,Es)

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
B=Vs(:,1); 
W=Vs(:,2); 

if(Es.JacMode==0)      % Model equations

    dB = Ps.lambda.*B.*W.*(1 + Ps.eta.*B).*(1 - B./Ps.kappa) - Ps.mu.*B; 
    dW = Ps.P - Ps.nu.*W.*(1-Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*W.*(1 + Ps.eta.*B);

    VsOut = [dB,dW];
else               % Jacobian of equations
    BdB = Ps.lambda.*W.* (1 + 2*B.*(Ps.eta - 1./Ps.kappa) - 3*Ps.eta./Ps.kappa.*B.^2) - Ps.mu;
    BdW = Ps.lambda.*B.*(1 + Ps.eta.*B).*(1 - B./Ps.kappa);
    WdB = - Ps.gamma.*W.*(1 + 2*B.*Ps.eta) + Ps.rho.*Ps.nu.*W./Ps.kappa;
    WdW = - Ps.nu.*(1-Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*(1 + Ps.eta.*B);

    % written in a large sparse matrix format 
    VsOut = ArrangeJacobian([BdB BdW;WdB WdW],Ps,Es);
end;





end