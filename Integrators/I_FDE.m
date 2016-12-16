function VsOut=I_FDE(Vs,Ps,Es,varargin)
% Integrator with Finite-Difference Euler scheme
% VsOut=I_FDE(Vs,Ps,Es)
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and non-local parts of the system

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

     % Setup the spatial matrix and auxiliary flags if not already done
    if(~isfield(Es,'SmUse'))
        [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);
    end;
    
    NonNeg = 0;		% If we know variables are positive, make sure they remain so	
    if((isfield(Es,'NonNeg')) & (Es.NonNeg))
        NonNeg = 1;
    end;
    
    totsteps = ceil(Es.TimeDst/Es.TsSize);
    % Go through each time step
    for ii=1:totsteps
        if(Es.SmUse)   % Integrate next time step
            %[ii size(Ps.LocFunc(Vs,Ps,Es)) size(Vs) size(Ps.SpaMat)]
            %size(reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.VarNum))
            VsNew = Vs + Es.TsSize*(Ps.LocFunc(Vs,Ps,Es) + reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.VarNum));
            %size(VsNew)
            if Es.SmUpdate
                Ps.SpaMat = Ps.SpaFunc(Vs,Ps,Es);  % Use this if the spatial matrix needs to be updated online
            end;
        else        % if we don't use SM (spatial matrix) than use the spatial function directly
            VsNew = Vs + Es.TsSize*(Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es));
        end;
    	if NonNeg  % consider deleting
        	VsNew = max(0,VsNew);
        end;
    	Vs = VsNew;
    end; 

    VsOut = Vs;



end
