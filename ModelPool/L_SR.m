function VsOut=L_SR(Vs,Ps,Es)
% simple model for mid-way between logistic and allee-effect dynamics - Local terms
% VsOut=L_SR(Vs,Ps,Es)
% Variables are: N(1) : dN/dt = r*N* (1-N/K)*(N/K)^gamma + d* D^2(N))
% Parameters are: r,K,gamma,d. (0.5,1,1,1)

% Initialization
N=Vs(:,1);

if(Es.JacMode==0)      % Model equations

    dN = Ps.r.*N.*(1 - N./Ps.K).*(N./Ps.K).^Ps.gamma;
    VsOut = dN;
else                % Jacobian of equations

    NdN = Ps.r.*(Ps.gamma+1 - (2+Ps.gamma)*N./Ps.K).*(N./Ps.K).^Ps.gamma;
    % written in a large sparse matrix format
    VsOut = spdiags(NdN,0,Ps.Nx*Ps.Ny,Ps.Nx*Ps.Ny);
end;

end
