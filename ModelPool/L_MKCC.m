function VsOut=L_MKCC(Vs,Ps,Es,varargin)
% Modified Klausmeier model with Carrying Capacity - Local terms
% VsOut=L_MKCC(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: B(1) and W(2). 
% Parameters are: a,m,K,D. (0.9,0.45,1,100)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initilization

B=Vs(:,1); 
W=Vs(:,2); 

if(Es.JacMode==0)      % Model equations

    dB = W.*B.*B.*(1-B./Ps.K) - Ps.m.*B;
    dW = Ps.a - W - W.*B.*B;

    VsOut = [dB,dW];
else                % Jacobian of equations
    BdB = W.*B.*(2-3*B./Ps.K) - Ps.m;
    BdW = B.*B;
    WdB = - 2*W.*B;
    WdW = - 1 - B.*B ;
    
    syslen = Ps.Nx * Ps.Ny;
    % written in a large sparse matrix format 
    VsOut = spdiags([BdB BdW WdB; WdW BdW WdB],[0 syslen -syslen],syslen*2,syslen*2);
    %VsOut = sparse([diag(BdB) diag(BdW) ; diag(WdB) diag(WdW)]);
end;


