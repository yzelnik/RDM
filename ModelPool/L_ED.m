function VsOut=L_ED(Vs,Ps,Es,varargin)
% Exponential decay - Local terms
% VsOut=L_ED(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Parameter is lambda for half-life

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
U=Vs(:,1); 

if(Es.JacMode==0)      % Model equations

    dU = -U./Ps.lambda;
    VsOut = [dU];
else               % Jacobian of equations

    % written in a large sparse matrix format 
    VsOut = ArrangeJacobian([(-1./Ps.lambda).*ones(size(U))],Ps,Es);
end;


end