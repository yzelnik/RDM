function VsOut=I_FDSIMP(Vs,Ps,Es,varargin)
% Integrator with Finite-Difference Semi-Implicit scheme
% VsOut=I_FDSIMP(Vs,Ps,Es)
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and spatial parts of the system
% Implicit scheme: U(T+t) = U(T) + t*L(U(T)) + t*D*U(T+t)
% where U represents the variables, T,t are the current time and time-step, 
% and L,S are the local and spatial functions (and 1 is the unity matrix)
% The scheme is rewritten  as: (1-t*D)*U(T+t) = (U(T)+t*L(U(T))) 
% So that U(T+t) is found by dividing (from the left) by the matrix (1-t*D)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'NonNeg',0);

% Setup the spatial matrix and auxiliary flags if not already done
if(~isfield(Es,'SpaMatUse') || Es.SpaMatUse<0)
	[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);
end;

if(~Es.SpaMatUse)
	error('Semi-Implicit integration is not possible without a matrix structure for the spatial part of the PDE.');
end;

% Calculate the matrix of (1-t*M) where 1 is unity matrix, t is time-step, and M is the spatial matrix.
tempmat = speye(size(Ps.SpaMat)) - Es.TsSize*Ps.SpaMat;
syslen = Ps.Nx * Ps.Ny;
totsteps = ceil(Es.TimeDst/Es.TsSize);
% Go through each time step
for ii=1:totsteps   
    % Explicit step for local (non-linear) part
    VsTemp = reshape( Vs + Es.TsSize*Ps.LocFunc(Vs,Ps,Es) ,syslen*Ps.VarNum,1);	
    % Implicit step for spatial (potentially linear) part
    Vs = reshape( tempmat \ VsTemp ,syslen,Ps.VarNum);    
    %VsChange = jac \ reshape(rhs,totlen,1);
    
	if Es.NonNeg  % make sure values are not negative, if relevant
    	Vs = max(0,Vs);
    end;
end; 

VsOut = Vs;

end
