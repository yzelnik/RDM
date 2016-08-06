
function VsOut=I_FDSIMP(Vs,Ps,Es,varargin)
% Integrator with Finite-Difference semi-Implicit scheme
% VsOut=I_FDSIMP(Vs,Ps,Es)
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and non-local parts of the system

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Setup the spatial matrix and auxiliary flags if not already done
if(~isfield(Es,'UseSM') )
	[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);
end;

if(~Es.UseSM)
	error('Semi-Implicit integration is not possible without a matrix strucutre for the spatial part of the PDE.');
end;

posflag = 0;		% If we know variables are positive, make sure they remain so	
if((isfield(Es,'posflag')) & (Es.posflag))
    posflag = 1;
end;

% Calculate the matrix of (1-t*M) where 1 is unity matrix, t is time-step, and M is the spatial matrix.
tempmat = speye(size(Ps.SpaMat)) - Es.Tstep*Ps.SpaMat;
syslen = Ps.Nx * Ps.Ny;
totsteps = ceil(Es.Tdest/Es.Tstep);
% Go through each time step
for ii=1:totsteps  
    % Explicit step for local (non-linear) part
    VsTemp = reshape( Vs + Es.Tstep*Ps.LocFunc(Vs,Ps,Es) ,syslen*Ps.Vnum,1);	
    % Implicit step for spatial (potentially linear) part
    Vs = reshape( tempmat \ VsTemp ,syslen,Ps.Vnum);    
    %VsChange = jac \ reshape(rhs,totlen,1);
    
	if posflag  % consider deleting
    	Vs = max(0,Vs);
    end;
    if Es.updateSM     % Use this if the spatial matrix needs to be updated online
        tempmat = speye(size(Ps.SpaMat)) - Es.Tstep*Ps.SpaFunc(Vs,Ps,Es);
    end;
end; 

VsOut = Vs;

end
