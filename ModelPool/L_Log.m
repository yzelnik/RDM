function VsOut=L_Log(Vs,Ps,Es,varargin)
% simple model for logistic dynamics - Local terms
% VsOut=L_Log(Vs,Ps,Es)
% Variables are: N(1) : dN/dt = r*N* (1-N/K) + d* D^2(N))
% Parameters are: r,K,d. (0.5,1,1)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
N=Vs(:,1); 

if(Es.JacMode==0)      % Model equations

    dN = Ps.r.*N.*(1-N./Ps.K).*(N./Ps.A-1);
    VsOut = dN;
else                % Jacobian of equations
    
    NdN = Ps.r .* ( -1+(2*N.*(Ps.A+Ps.K)-3*N.^2)./(Ps.A.*Ps.K) );
    % written in a large sparse matrix format 
    VsOut = spdiags(NdN,0,Ps.VarNum,Ps.VarNum);
end;

end
