function VsOut=I_NoiseEM(Vs,Ps,Es,varargin)
% Integrator with noise (SDE) using the Euler-Maruyama scheme
% VsOut=I_NoiseEM(Vs,Ps,Es)
% Ps.NoiseType should be either a number, a function or a cell array:
% - number(s):  [a b] --> noise(x) = a*x.^b  (where x is the variable)
% - function:   a function handle, assuming no extra input is needed
% - cell-array: {func,prms} where func is some function, and prms its parameters
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and non-local parts of the system

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'NonNeg',0);

% Setup the spatial matrix and auxiliary flags if not already done
if(~isfield(Es,'SpaMatUse') || Es.SpaMatUse<0)
    [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);
end;

% setup noise funtion and parameters
[noisefunc,noiseparm]=SetupNoise(Ps,Es);

% how many steps of integration?
totsteps = ceil(Es.TimeDst/Es.TsSize);
    
% Go through each time step
for ii=1:totsteps
	for jj=1:Ps.VarNum % setup noise
        tmpnoiseform(:,jj) = noisefunc(Vs(:,jj),noiseparm(:,jj));
	end;
	noise = sqrt(Es.TsSize).*randn(size(Vs)).*tmpnoiseform;
	if(Es.SpaMatUse)   % Integrate next time step
        VsNew = Vs + Es.TsSize*(Ps.LocFunc(Vs,Ps,Es) + reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.VarNum)) + noise;
	else        % if we don't use SM (spatial matrix) than use the spatial function directly
        VsNew = Vs + Es.TsSize*(Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es)) + noise;
	end;
	if Es.NonNeg % make sure values are not negative, if relevant
        VsNew = max(0,VsNew);
	end;
	Vs = VsNew;
end; 
VsOut = Vs;

end

