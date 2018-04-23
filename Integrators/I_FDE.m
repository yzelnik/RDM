function VsOut=I_FDE(Vs,Ps,Es,varargin)
% Integrator with Finite-Difference Euler scheme
% VsOut=I_FDE(Vs,Ps,Es)
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and non-local parts of the system

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'NonNeg',0);

% Setup the spatial matrix and auxiliary flags if not already done
if(~isfield(Es,'SpaMatUse') || Es.SpaMatUse<0)
	[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);
end;

totsteps = ceil(Es.TimeDst/Es.TsSize);

% Go through each time step
for ii=1:totsteps
	if(Es.SpaMatUse)   % Integrate next time step
        VsNew = Vs + Es.TsSize*(Ps.LocFunc(Vs,Ps,Es) + reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.VarNum));
	else  % if we don't use SM (spatial matrix) than use the spatial function directly
        VsNew = Vs + Es.TsSize*(Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es));
	end;
	if Es.NonNeg  % make sure values are not negative, if relevant
        VsNew = max(0,VsNew);
	end;        
    Vs = VsNew;
end; 
VsOut = Vs;

end


