function VsOut=M_InitPerSt(Vs,Ps,Es,varargin)
% Initiate a periodic state/perturbation using Vs
% VsOut=M_InitPerState(Vs,Ps,Es)
% For 1D - stripes; For 2D - either stripes of hexagons;
% Assumes: Es.InitParms = [NoPer, PerType]
% NoPer is simply the period number, while PerType: 0=stripes, 1=spots, (-1)=holes;
% If only one uniform state (in Vs) is supplied, a perturbation of the pattern is made
% Otherwise, the patterned is formed from the first 2 uniform states - Vs(:,:,1:2)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'InitParms'))
	Es.InitParms = [4 1 0];
end;

NoPer = Es.InitParms(1);	% Period number
PerType = Es.InitParms(2);	% Pattern type: 0=stripes, 1=spots, (-1)=holes

% prepare x and k
x = linspace(0,Ps.Lx,Ps.Nx)';
lambda = Ps.Lx/NoPer;
k = 2*pi/lambda;

% check if this is a 1D system
if((Ps.Nx==1) || (Ps.Ny==1))
	u = (1+cos(k.*x))/2;
	u = repmat(u,1,Ps.Vnum);
else  % Assuming this is a 2D system
	%[y,x] = meshgrid(linspace(0,Ps.Ly,Ps.Ny),x);
	y = linspace(0,Ps.Ly,Ps.Ny);
	[x,y] = meshgrid(x,y);
	
	% Build initial mask
	u = (cos(k.*x) + abs(PerType) * (cos(-k/2.*x+sqrt(3)/2*k.*y) + cos(-k/2.*x-sqrt(3)/2*k.*y)));
	% Normalize to: 0<u<1
	mn = min(u(:)); 
	mx = max(u(:));
	u = (u-mn)/(mx-mn);
	if(PerType<0)	% Flip u (to get holes, etc)
		u = 1-u;
	end;
	%imagesc(u)
	u = repmat(reshape(u',length(u(:)),1),1,Ps.Vnum);
	
end

if(size(Vs,1)<(Ps.Nx*Ps.Ny))
	if(size(Vs,1)>2)	% Cut the Vs short if it's too long
		Vs=Vs(1:2,:);
	end;
	% Prepare Uniform states if necessary
	Vs = M_InitUnfSt(Vs,Ps,Es);
end;

if (size(Vs,3)==1)	% If only one set of Vs is specified, make a perturbation
	VsOut = Vs + (u*2-1)*Es.STsmall; 
else			% Use first 2 sets of Vs to build an initial state
	VsOut = u .* Vs(:,:,1) + (1-u) .* Vs(:,:,2);
end;

% If we know variables are positive, make sure they remain so	
if((isfield(Es,'posflag')) & (Es.posflag))
	VsOut = max(0,VsOut);
end;

end
