function VsOut=IFD_RK4(Vs,Ps,Es,varargin)
% Integrator with Finite-Difference - using Runge-Kutta 4th order
% VsOut=IFD_RK4(Vs,Ps,Es)
% LocFunc and SpaFunc are the functions that detail the local and non-local parts of the system

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

SpaMat = Ps.SpaFunc(Vs,Ps,Es);
% Is the result of the SpaFunc a SM, or just a "state" vector?
if(size(SpaMat,1)==size(SpaMat,2))
	UseSM=1;
else
	UseSM=0;
end;

posflag=0;		% If we know variables are positive, make sure they remain so	
if((isfield(Es,'posflag')) & (Es.posflag))
	posflag = 1;
end;
updateSM=0;		% Check if the Spatial matrix needs to be updated (due to Non-Linear Derivatives)
if((isfield(Ps,'NLD')) & (Ps.NLD))
	updateSM = 1;
end;

% Go through each time step
for ii=1:ceil(Es.Tdest/Es.Tstep)
	% Integrate next time step
	if(UseSM)
		F1 = (Ps.LocFunc(Vs,Ps,Es) + reshape(SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.Vnum));
		VF1= Vs + F1*Es.Tstep/2;
		F2 = (Ps.LocFunc(VF1,Ps,Es) + reshape(SpaMat*VF1(:),Ps.Nx*Ps.Ny,Ps.Vnum));
		VF2= Vs + F2*Es.Tstep/2;
		F3 = (Ps.LocFunc(VF2,Ps,Es) + reshape(SpaMat*VF2(:),Ps.Nx*Ps.Ny,Ps.Vnum));
		VF3= Vs + F3*Es.Tstep;
		F4 = (Ps.LocFunc(VF3,Ps,Es) + reshape(SpaMat*VF3(:),Ps.Nx*Ps.Ny,Ps.Vnum));
		VsNew = Vs + Es.Tstep*(F1 + 2*F2 + 2*F3 + F4)/6;
		if updateSM
			SpaMat = Ps.SpaFunc(Vs,Ps,Es);
		end;
	else
		VsNew = Vs + Es.Tstep*(Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es));
	end;
	if posflag
		VsNew=max(0,VsNew);
	end;
	
	Vs = VsNew;
end; 

VsOut = Vs;

end
