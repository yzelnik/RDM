function VsOut = M_InitCorrRnd(Vs,Ps,Es,varargin)
% Form a new correlated-random state around a uniform value (Given by M_InitUnfSt)
% VsOut = M_InitCorrRnd(Vs,Ps,Es,varargin)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'NonNeg',0,'VarInd',1);

% Prepare Uniform states if necessary
VsOut = M_InitUnfSt(Vs,Ps,Es);

% what is the nx/ny dimension? (how many pixels)
if(Ps.Ny==1)
    syslen=Ps.Nx;
elseif(Ps.Nx==1)
    syslen=Ps.Ny;
else
    syslen=(Ps.Nx+Ps.Ny)/2;
end;
% what is the correation size to try for?
if(Es.CorrSize<1)
    Es.CorrSize=Es.CorrSize*syslen;
end;
Es.CorrSize = syslen/Es.CorrSize;
randmask = rand(Ps.Nx,Ps.Ny);

% go to frequency domain, and take out high frequencies
tmpfft = fftshift(fft2(randmask));
[xx,yy]=meshgrid((0.5:Ps.Ny)-Ps.Ny/2,(0.5:Ps.Nx)-Ps.Nx/2);
circ = sqrt(xx.^2+yy.^2)<Es.CorrSize;

tmpout = abs(ifft2(fftshift(circ.*tmpfft)));
randmask = tmpout(:)-min(tmpout(:));
randmask = randmask/max(randmask);

% Add a perturbation around the uniform state
VsOut(:,Es.VarInd) = randmask*2*VsOut(1,Es.VarInd);

if(Es.NonNeg)  % make sure values are not negative, if relevant
	VsOut = max(0,VsOut);
end;

end
