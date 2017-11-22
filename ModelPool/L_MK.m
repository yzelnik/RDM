function VsOut=L_MK(Vs,Ps,Es,varargin)
% Modified Klausmeier model - Local terms
% VsOut=L_MK(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1) and W(2). 
% Parameters are: a,m,D. (0.9,0.45,100)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initilization

B=Vs(:,1); 
W=Vs(:,2); 

if(Es.JacMode==0)      % Model equations

    dB = W.*B.*B - Ps.m.*B;
    dW = Ps.a - W - W.*B.*B;

    VsOut = [dB,dW];
else                % Jacobian of equations
    BdB = 2*W.*B - Ps.m;
    BdW = B.*B;
    WdB = - 2*W.*B;
    WdW = - 1 - B.*B ;
    
    %syslen = Ps.Nx * Ps.Ny;
    % written in a large sparse matrix format 
    VsOut = ArrangeJacobian([BdB BdW;WdB WdW],Ps,Es);
end;


