function VsOut=IFD(Vs,Ps,Es,varargin)
% Integrator with Finite-Difference Euler scheme
% VsOut=IFD(Vs,Ps,Es)
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and non-local parts of the system

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

%if(isfield(Ps,'SpaMat'))
%    SpaMat=Ps.SpaMat;
%else
    SpaMat = Ps.SpaFunc(Vs,Ps,Es);
%end;
% Is the result of the SpaFunc a SM, or just a "state" vector?
if(size(SpaMat,1)==size(SpaMat,2))
	UseSM = 1;
else
	UseSM = 0;
end;

posflag = 0;		% If we know variables are positive, make sure they remain so	
if((isfield(Es,'posflag')) & (Es.posflag))
	posflag = 1;
end;
updateSM = 0;		% Check if the Spatial matrix needs to be updated (due to Non-Linear Derivatives)
if((isfield(Ps,'NLD')) & (Ps.NLD))
	updateSM = 1;
end;

% Go through each time step
for ii=1:ceil(Es.Tdest/Es.Tstep)
    %disp([ii ceil(Es.Tdest/Es.Tstep)]);
    
	% Integrate next time step
	if(UseSM)
		VsNew = Vs + Es.Tstep*(Ps.LocFunc(Vs,Ps,Es) + reshape(SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.Vnum));
		if updateSM
			SpaMat = Ps.SpaFunc(Vs,Ps,Es);
		end;
    else
        
		VsNew = Vs + Es.Tstep*(Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es));
        
	end;
% consider deleting
	if posflag
		VsNew = max(0,VsNew);
	end;
	
	Vs = VsNew;
    
end; 

VsOut = Vs;

end
