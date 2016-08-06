function VsOut=InitPerState(Vs,Ps,Es,NoPer,PerType,varargin)
% Initiate a periodic perturbation on top of Vs
% For 1D - stripes; For 2D - either stripes of hexagons;
% NoPer = number of periods; 
% PerType: 'str'=stripes, 'spo'=hexagonal spots, 'hol'=hexagonal holes;
% VsOut=InitPerState(Vs,Ps,Es,NoPer,PerType)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if(nargin<4)
	NoPer=4;
end;
if(nargin<5)
	PerType='spo';
end;

% prepare x and k
x = linspace(0,Ps.Lx,Ps.Nx)';
lambda = Ps.Lx/NoPer;
k = 2*pi/lambda;

% check if this is a 1D system
if((Ps.Nx==1) || (Ps.Ny==1))
    u = cos(k.*x);

else  % Assuming this is a 2D system
	%[y,x] = meshgrid(linspace(0,Ps.Ly,Ps.Ny),x);
	[x,y] = meshgrid(x,linspace(0,Ps.Ly,Ps.Ny));
	%size(x)
	%size(y)
	if strcmp(PerType,'str')        
		u = cos(k.*x);
	elseif strcmp(PerType,'spo')        
		u = cos(k.*x)+cos(-k/2.*x+sqrt(3)/2*k.*y)+cos(-k/2.*x-sqrt(3)/2*k.*y);	
	elseif strcmp(PerType,'hol')        
		u = cos(k.*x)+cos(-k/2.*x+sqrt(3)/2*k.*y)+cos(-k/2.*x-sqrt(3)/2*k.*y);	
		u = 1-u;
	end;
	u = reshape(u,length(u(:)),1);
	
end

% prepare output to avoid negative values
%u = u + abs(min(min(u)));
u = Es.STsmall.*u;
Vper = repmat(u,1,Ps.Vnum);
%size(Vs)
%size(Vper)
VsOut = Vs + Vper;

end
