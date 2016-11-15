function [VsOut,residual,ind]=NewtonLoop(Vs,Ps,Es,varargin)
% Repeat a Newton loop until threshold is reached
% [VsOut,ind,residual]=NewtonLoop(Vs,Ps,Es)
% Ps.LocFunc & Ps.SpaFunc are the functions for the local & non-local parts
% The threshold is set by Es.SsThresh (or 1e-10 by default)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Setup the spatial matrix and auxiliary flags
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es,varargin{:});

Es=InsertDefaultValues(Es,'MaxNewtLoop',20,'NonNeg',0,'SsThresh',1e-10,'JacNum',0,'NewtonRate',1,'InverseMatrix',0);

syslen = Ps.Nx * Ps.Ny;
totlen = syslen * Ps.VarNum;

ind = 1;
residual = inf;
while (ind<Es.MaxNewtLoop) && (residual>Es.SsThresh)
    % Get the right-hand-side of the PDE (no time derivatives)
    rhs=RightHandSide(Vs,Ps,Es);  
    % Calculate the residuals, and normalize
    residual = (sqrt(sum(rhs(:).^2))/totlen);
    
    jac=CalculateJacobian(Vs,Ps,Es);
    
    % Use the inverse of the jac matrix on the rhs vector to update the state vector
    if(Es.InverseMatrix) % find inverse of matrix?
        disp(999)
        invjac = inv(jac);
        VsChange = invjac * reshape(rhs,totlen,1); 
    else  % or use the faster direct matlab division operator
        VsChange = jac \ reshape(rhs,totlen,1); 
    end;
    
    Vs = Vs - Es.NewtonRate*reshape(VsChange,syslen,Ps.VarNum);

    if Es.NonNeg  % consider deleting
        	Vs = max(0,Vs);
    end;
    ind = ind+1;

end;

VsOut = Vs;
%disp([ind log10(residual)])% delete this line

end
