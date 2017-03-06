function VsOut=I_NoiseEM(Vs,Ps,Es,varargin)
% Integrator with noise (SDE) using the Euler-Maruyama scheme
% VsOut=I_NoiseEM(Vs,Ps,Es)
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and non-local parts of the system

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'NonNeg',0);

% Setup the spatial matrix and auxiliary flags if not already done
if(~isfield(Es,'SmUse'))
    [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);
end;

% setup noise funtion and parameters
[noisefunc,noiseparm]=SetupNoise(Ps,Es);

% how many steps of integration?
totsteps = ceil(Es.TimeDst/Es.TsSize);
    
    %disp([totsteps Es.TsSize Es.TimeDst])
    % Go through each time step
    for ii=1:totsteps
        for jj=1:Ps.VarNum % setup noise
            tmpnoiseform(:,jj) = noisefunc(Vs(:,jj),noiseparm(:,jj));
        end;
        noise = sqrt(Es.TsSize).*randn(size(Vs)).*tmpnoiseform;
        
        %noise = sqrt(Es.TsSize).*randn(size(Vs)).*Vs/2;
        %mean(abs(noise(:)))
        if(Es.SmUse)   % Integrate next time step
            %[ii size(Ps.LocFunc(Vs,Ps,Es)) size(Vs) size(Ps.SpaMat)]
            %size(reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.VarNum))
            VsNew = Vs + Es.TsSize*(Ps.LocFunc(Vs,Ps,Es) + reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.VarNum)) + noise;
            %size(VsNew)
            if Es.SmUpdate
                Ps.SpaMat = Ps.SpaFunc(Vs,Ps,Es);  % Use this if the spatial matrix needs to be updated online
            end;
        else        % if we don't use SM (spatial matrix) than use the spatial function directly
            VsNew = Vs + Es.TsSize*(Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es)) + noise;
        end;
    	if Es.NonNeg  
        	VsNew = max(0,VsNew);
        end;
    	Vs = VsNew;
    end; 

    VsOut = Vs;

end
