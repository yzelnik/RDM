function VsOut=I_NoiseDFM(Vs,Ps,Es,varargin)
% Integrator with noise (SDE) using the derivative free Milstein scheme
% VsOut=I_NoiseDFM(Vs,Ps,Es)
% Ps.NoiseType should be either a number, a function or a cell array:
% - number(s):  [a b] --> noise(x) = a*x.^b  (where x is the variable, b=0 by default)
% - function:   a function handle, assuming no extra input is needed
% - cell-array: {func,prms} where func is some function, and prms its parameters
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and non-local parts of the system
% X_next = X + f(X)*dt + g(X)*W + (g(Y)-g(X))*(W^2)/(2*sqrt(dt)) 
% where f,g are the rhs without noise and the noise function, and X are the varibles
% Y = X + f(X)*dt + g(x)*sqrt(dt)
% W is the noise term, ~ sqrt(dt) * Gauss(0,1)

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
	% Get right-hand-side (without noise)
	if(Es.SpaMatUse)   % Integrate next time step   
        fx = (Ps.LocFunc(Vs,Ps,Es) + reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.VarNum));
	else  % if we don't use SM (spatial matrix) than use the spatial function directly
        fx = (Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es)) ;
	end;

    for jj=1:Ps.VarNum % calculate noise for this step
        gx(:,jj) = noisefunc(Vs(:,jj),noiseparm(:,jj));
        nextvar = Vs(:,jj) + fx(:,jj)*Es.TsSize + gx(:,jj)*sqrt(Es.TsSize);
        gnext(:,jj) = noisefunc(nextvar,noiseparm(:,jj));
	end;
	Wnoise = sqrt(Es.TsSize)*randn(size(Vs));
   
	% Next time step
	VsNew  = Vs + fx*Es.TsSize + gx.*Wnoise + (gnext-gx).*(Wnoise.^2)/(2*sqrt(Es.TsSize));
        
	if Es.NonNeg % make sure values are not negative, if relevant
        VsNew = max(0,VsNew);
	end;
	Vs = VsNew;
end; 
VsOut = Vs;

end

