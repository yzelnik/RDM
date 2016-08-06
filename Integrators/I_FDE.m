function VsOut=I_FDE(Vs,Ps,Es,varargin)
% Integrator with Finite-Difference Euler scheme
% VsOut=I_FDE(Vs,Ps,Es)
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and non-local parts of the system

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'fmod'))
   Es.fmod=0;
end;

if(Es.fmod<0)   % Setup variables for integration in the future
    
    Ps = Ps.SpaFunc(Vs,Ps,Es);    % Get spatial matrix for future integration
    VsOut = Ps; % This is a "misuse" of the name, but we just want to return the "new" Ps struct
    
else            % Normal run
     % Setup the spatial matrix and auxiliary flags if not already done
    if(~isfield(Es,'UseSM'))
        [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);
    end;
    
    posflag = 0;		% If we know variables are positive, make sure they remain so	
    if((isfield(Es,'posflag')) & (Es.posflag))
        posflag = 1;
    end;
    
    totsteps = ceil(Es.Tdest/Es.Tstep);
    % Go through each time step
    for ii=1:totsteps
        if(Es.UseSM)   % Integrate next time step
            VsNew = Vs + Es.Tstep*(Ps.LocFunc(Vs,Ps,Es) + reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.Vnum));
            if Es.updateSM
                Ps.SpaMat = Ps.SpaFunc(Vs,Ps,Es);  % Use this if the spatial matrix needs to be updated online
            end;
        else        % if we don't use SM (spatial matrix) than use the spatial function directly
            VsNew = Vs + Es.Tstep*(Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es));
        end;
    	if posflag  % consider deleting
        	VsNew = max(0,VsNew);
        end;
    	Vs = VsNew;
    end; 

    VsOut = Vs;

end;        % End of normal run

end
