function eigvals=T_SpatialDynamicsEigenvalues(Vs,Ps,Es,varargin)
% Use the Jacobian to find the eigenvalues for the spatial dynamics of uniform states
% eigs=T_SpatialDynamicsEigenvalues
% A reaction-diffusion setup is currently assumed, so number of equations
% is double that of the number of variables
% Uses the Ps.LocFunc to construct a Jacobian

Es.InitActive = 1;	% If only the ODE state is given, don't look for more 
% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es.JacMode = 1;	% Request a jacobian

if(size(Vs,3)>1)    % We expect a set of uniform states, where the first axis runs through them.
    UnifStates = squeeze(Vs(1,:,:))';
else
    if((Ps.Nx*Ps.Ny)==size(Vs,1))
        UnifStates = Vs(1,:); % It is assumed just one state was given, and should be treated as uniform
    else
        UnifStates = Vs;
    end;
end;
Ps.Nx=1; 
Ps.Ny=1;

for ii=1:size(UnifStates,1)
    locjac = full(Ps.LocFunc(UnifStates(ii,:),Ps,Es));
    spadynmat = [zeros(Ps.VarNum) eye(Ps.VarNum) ; repmat(-1./Ps.Ds(:),1,Ps.VarNum).*locjac zeros(Ps.VarNum) ];
    eigvals(ii,:) = transpose(eigs(spadynmat));
end;

end

 

