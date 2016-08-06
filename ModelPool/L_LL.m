function VsOut=L_LL(Vs,Ps,Es,varargin)
% Lejeune-Lefever model - Local terms
% VsOut=L_LL(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1)
% Parameters are: mu,Lambda,L. (1.02,1.2,0.2)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if(~isfield(Es,'fmod'))
   Es.fmod=0;
end;

% Initialization

B = Vs(:,1); 

if(Es.fmod==0)      % Model equations

    dB = (1 - Ps.mu).*B + (Ps.Lambda - 1).*B.^2 - B.^3;
    VsOut = dB;
else                % Jacobian of equations

    BdB = 1 - Ps.mu + 2*(Ps.Lambda - 1).*B - 3*B.^2;
    %VsOut = sparse(diag(BdB));  % written in a large sparse matrix format 
    VsOut = spdiags(BdB,0,length(BdB),length(BdB)); % written in a large sparse matrix format 
end;

end
