function VsOut=L_SR(Vs,Ps,Es,varargin)
% Slow recovery model, with "deep hole" at N=K, but "flat surface" outside - Local terms
% VsOut=L_SR(Vs,Ps,Es)
% Variables are: N(1) : dN/dt = -r*alpha*sinh(alpha*N-alpha*K)/cosh(alpha*N-alpha*K).^2 + d* D^2(N))
% "energy" functional looks like 1-r/cosh(alpha(N-K)) so that there is a hole at N=K,
% where outside of it the energy-landscape becomes flat very quickly (for large alpha)
% Parameters are: r,K,alpha,d. (1,1,10,1)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
N=Vs(:,1); 

if(Es.JacMode==0)      % Model equations
    dN = -Ps.r*Ps.alpha.*sinh(Ps.alpha.*N-Ps.alpha.*Ps.K)./(cosh(Ps.alpha.*N-Ps.alpha.*Ps.K).^2);
    VsOut = dN;
else                % Jacobian of equations
    
    NdN = 0; 
    %warning('jacobian not implemented');
    % written in a large sparse matrix format 
    VsOut = spdiags(NdN,0,1,1);
end;

end
