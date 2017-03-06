function VsOut=L_GKGS(Vs,Ps,Es,varargin)
% Generalized Klausmeier-Gray-Scott model - Local terms
% VsOut=L_GKGS(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1) and W(2), Parameters are: a,m,c,D. (0.9,0.45,1,100)
% Equations: dB/dt = W*B^2 - m*B + B'',  dW/dt = a - c*W - W*B^2 + D*W''

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initilization
B=Vs(:,1); 
W=Vs(:,2); 

if(Es.JacMode==0)      % Model equations
    
    dB = W.*B.^2 - Ps.m.*B;
    dW = Ps.a - Ps.c*W - W.*B.^2;
    
    VsOut = [dB,dW];
else    % Jacobian of equations
    
    BdB = 2*W.*B - Ps.m;
    BdW = B.^2;
    WdB = - 2*W.*B;
    WdW = - Ps.c - B.^2 ;
    
    % written in a large sparse matrix format 
    VsOut = ArrangeJacobian([BdB BdW;WdB WdW],Ps,Es);
end;


