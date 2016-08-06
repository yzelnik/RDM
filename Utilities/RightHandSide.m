function rhs=RightHandSide(Vs,Ps,Es,varargin)
% Return right-hand-side of equation (both local and spatial parts)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Setup the spatial matrix and auxiliary flags
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);		

if(Es.UseSM)   % Integrate next time step 
    if Es.updateSM
        Ps.SpaMat = Ps.SpaFunc(Vs,Ps,Es);  % Use this if the spatial matrix needs to be updated online        
    end;
    rhs = (Ps.LocFunc(Vs,Ps,Es) + reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.Vnum));     
else        % if we don't use SM (spatial matrix) than use the spatial function directly
    rhs = (Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es));
end;

end
