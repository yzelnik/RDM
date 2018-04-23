function VsOut=I_FDRK4(Vs,Ps,Es,varargin)
% Integrator with Finite-Difference Runge-Kutta 4th order scheme
% VsOut=I_FDRK4(Vs,Ps,Es)
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
        F1 = (Ps.LocFunc(Vs,Ps,Es) + reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.VarNum));
        VF1= Vs + F1*Es.TsSize/2;
        F2 = (Ps.LocFunc(VF1,Ps,Es) + reshape(Ps.SpaMat*VF1(:),Ps.Nx*Ps.Ny,Ps.VarNum));
        VF2= Vs + F2*Es.TsSize/2;
        F3 = (Ps.LocFunc(VF2,Ps,Es) + reshape(Ps.SpaMat*VF2(:),Ps.Nx*Ps.Ny,Ps.VarNum));
        VF3= Vs + F3*Es.TsSize;
        F4 = (Ps.LocFunc(VF3,Ps,Es) + reshape(Ps.SpaMat*VF3(:),Ps.Nx*Ps.Ny,Ps.VarNum));
        VsNew = Vs + Es.TsSize*(F1 + 2*F2 + 2*F3 + F4)/6;
	else        % if we don't use SM (spatial matrix) than use the spatial function directly
        F1 = (Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es));
        VF1= Vs + F1*Es.TsSize/2;
        F2 = (Ps.LocFunc(VF1,Ps,Es) + Ps.SpaFunc(VF1,Ps,Es));
        VF2= Vs + F2*Es.TsSize/2;
        F3 = (Ps.LocFunc(VF2,Ps,Es) + Ps.SpaFunc(VF2,Ps,Es));
        VF3= Vs + F3*Es.TsSize;
        F4 = (Ps.LocFunc(VF3,Ps,Es) + Ps.SpaFunc(VF3,Ps,Es));
        VsNew = Vs + Es.TsSize*(F1 + 2*F2 + 2*F3 + F4)/6;
	end;
    	
    if Es.NonNeg  % make sure values are not negative, if relevant
        VsNew = max(0,VsNew);
	end;
	Vs = VsNew;
 end; 
VsOut = Vs;

end

