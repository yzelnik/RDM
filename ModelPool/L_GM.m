function VsOut=L_GM(Vs,Ps,Es,varargin)
% Gierer-Meinhardt model with Saturation - Local terms
% VsOut=L_GM(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: a(1) and h(2), Parameters are: rho,mu,del,eps,D. (0.1,0.5,0.1,0.01,10)
% Equations: da/dt = rho + a^2/(h*(1+del*a^2)) - mu*a + a'',  dh/dt = eps + a^2 - h + D*U''
% Formulation based on paper: "Analytical Treatment of Pattern Formation
% in the Gierer-Meinhardt Model of Morphogenesis" 

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if(~isfield(Es,'fmod'))
   Es.fmod=0;
end;

% Initilization
a=Vs(:,1); 
h=Vs(:,2); 

if(Es.fmod==0)      % Model equations
    
    da = Ps.rho + a.^2./(h.*(1+Ps.del.*a.^2)) - Ps.mu.*a;
    dh = Ps.eps + a.^2 - h;
    
    VsOut = [da,dh];
else    % Jacobian of equations
    
    ada = + 2*a./(h.*(1+Ps.del.*a.^2).^2) - Ps.mu;
    adh = - a.^2./h.^2;
    hda = + 2*a;
    hdh = - 1;
    
    % written in a large sparse matrix format 
    VsOut = sparse([diag(ada) diag(adh) ; diag(hda) diag(hdh)]);
end;


