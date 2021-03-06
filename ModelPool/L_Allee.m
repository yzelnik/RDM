function VsOut=L_Allee(Vs,Ps,Es)
% simple model for the allee effect - Local terms
% VsOut=L_Allee(Vs,Ps,Es)
% Variables are: N(1) : dN/dt = r*N* (1-N/K)*(N/A-1) + d* D^2(N))
% Parameters are: r,K,A,d. (0.5,1,0.2,1)

% Initialization
N=Vs(:,1); 

if(Es.JacMode==0)      % Model equations

    dN = Ps.r.*N.*(1-N./Ps.K).*(N./Ps.A-1);%.*Ps.K.*Ps.A.*(Ps.K-Ps.A).^-2;
    VsOut = dN;
else                % Jacobian of equations
    
    NdN = Ps.r .* ( -1+(2*N.*(Ps.A+Ps.K)-3*N.^2)./(Ps.A.*Ps.K) );
    % written in a large sparse matrix format 
    VsOut = spdiags(NdN,0,Ps.Nx*Ps.Ny,Ps.Nx*Ps.Ny);
end;

end
