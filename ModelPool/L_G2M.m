function VsOut=L_G2M(Vs,Ps,Es,varargin)
% General 2 mechanism model - Local terms
% VsOut=L_G2M(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: U(1), V(2) and W(3), Parameters are: a,m,c,D. (0.9,0.45,1,100)
% Equations: dB/dt = W*B^2 - m*B + B'',  dW/dt = a - c*W - W*B^2 + D*W''
%Ps_G2M=struct('LocFunc',@L_G2M,'SpaFunc',@S_RD,'p',5,'a',1.5,'b',2,'c',1,'d',3,'e',4,'f',1,'g',1,'h',5,'i',3,'j',2,'Ds',[1 100 10000 100],'Vnum',4,'Lx',100,'Ly',1,'Nx',1000,'Ny',1);

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Initilization

%U=Vs(:,1); 
%V=Vs(:,2); 
%W=Vs(:,3); 

B=Vs(:,1); 
W=Vs(:,2); 
H=Vs(:,3); 
C=Vs(:,4); 


if(Es.Jcob==0)
% Model equations

%dU = Ps.a.* W./V.*U.^2 - Ps.b.*U;
%dV = Ps.c.* U.^2 - Ps.d.*V - Ps.g.*W ;
%dW = Ps.e  - Ps.f.*W.*U.^2 + Ps.g.*V ;

dB = Ps.a.* W.*B.*(1+Ps.e.*B) - Ps.b.*B;
dW = Ps.h.* H.*(1+Ps.y.*B) -Ps.c.*W - Ps.d.* W.*B.*(1+Ps.e.*B);
dH = Ps.p - Ps.h.* H.*(1+Ps.y.*B) - Ps.i.* H.*C.*(1+Ps.j.*C);
dC = Ps.f.* H.*C.*(1+Ps.j.*C) - Ps.g.*C;

VsOut = [dB,dW,dH,dC];


%VsOut = [dU,dV,dW];

else
% Jacobian of equations
UdU = 2*W.*U - Ps.m;
UdW = U.^2;
WdU = - 2*W.*U;
WdW = - Ps.c - U.^2 ;

% written in a large sparse matrix format 
VsOut = sparse([diag(UdU) diag(UdW) ; diag(WdU) diag(WdW)]);
end;


