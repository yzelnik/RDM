function [K,Lambda,Phi] = T_CalcWL(Vs,Ps,Es,varargin)
% Calculate wavenumber (k) and wavelength (Lambda) of a given state
% [K,Lambda,Phi] = T_CalcWL(Vs,Ps,Es)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if (~isfield(Es,'VarInd'))
    Es.VarInd = 1;
end;

if((Ps.Nx==1) || (Ps.Ny==1))
	dim=1;
else
	dim=2;
end;

vmax = max(Vs(:,Es.VarInd(1),1));
vmin = min(Vs(:,Es.VarInd(1),1));
if((vmax-vmin)<Es.StSmall)
	K=0;
	Lambda=0;
	Phi=zeros(1,dim);
else

if(dim==1)
	Phi = 0;		% For consistency with 2D
	U1 = fft(Vs(:,Es.VarInd(1),1));                                 % Fourier transform of Vs
	
    K = FindProminentMax(abs(U1(2:Ps.Nx/2)))*(2*pi)/Ps.Lx;
    %plot(Uabs)
    Lambda = 2*pi/K;                                       % wavelength
else	% dim==2
	V2 = reshape(Vs(:,Es.VarInd(1),1),Ps.Nx,Ps.Ny);		% Treating it as a 2D matrix
	U2 = fft2(V2);					% 2D FFT
	U2(1,1) = 0;

    U2 = abs(fftshift(fft2(V2)));   % Apply Fourier Transform, look at absolute values only
    [yy,xx]=meshgrid((1:Ps.Ny)-Ps.Ny/2-1,(1:Ps.Nx)-Ps.Nx/2-1);  
    ratio = Ps.Nx/Ps.Ny;    % assuming the resolution per pixel is the same
    dist = round(sqrt(xx.^2+(ratio*yy).^2));    % distance of each point from center
    for ii=1:ceil(min(Ps.Nx,Ps.Ny)/2-1)     % Divided by 3 (and not 2) as resolution is problematic in any case
        inds = find(dist==ii);              % Pixel locations for this distance
        rho(ii)=sum(U2(inds))/length(inds); % Averaging values at these locations
    end;
    
    maxloc = FindProminentMax(rho);
    K = maxloc*(2*pi)/Ps.Lx;
    Lambda = 2*pi/K; 
    
    inds = find((dist>maxloc-1).*(dist<maxloc+1));       % looking at points in the most dominant distance
    vals = U2(inds);                        % values at these points
    [xi,yi]=ind2sub([Ps.Nx,Ps.Ny],inds);    % relative locations (in x,y) of these points

    degs = atan((yi-Ps.Ny/2-1)./(xi-Ps.Nx/2-1));    % arctan of each of these points
    degs (degs>pi/2-eps) = -pi/2;                   % fixing a tan->+-pi/2 problem
    Phi  = -(degs'*vals)/sum(vals);          % averaging arctan's
    %[degs xi-Ps.Nx/2-1 yi-Ps.Ny/2-1 vals/100]
end;

end;


end


%%%%%%%%% AUX FUNCTION %%%%%%%%%%%%%%%%
function maxloc = FindProminentMax(vect)
% This function takes a series (signal), and tries to find prominent maximum
% Mainly, one tries to both ignore things right at the begining, and to estimate
% not a single integer location, but some effective average of its surroundings
vect=vect(:);
if(isempty(vect))
    maxloc=NaN;
else
    %plot(vect)
    reg = 2; % how much to look around max to estimate the effective location
    vectsm=smooth(vect); % Smoothing out the signal

    curv = diff(diff(vectsm))<0;    % negative curvature
    
    curv = ([0;0;curv] + [0;curv;0] + [curv;0;0])>1; % look at regions of negative curvature - where max should be
    
    [~,loc]=max((vect(:)+vectsm(:)).*curv); % find max value in negative curvature region
    inds = max(1,loc-reg):min(loc+reg,length(vect)); % locations used to average max
    maxloc = loc+((inds-loc)*vect(inds))/sum(vect(inds)); % final assessment

end;
end