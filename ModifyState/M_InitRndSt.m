function VsOut = M_InitRndSt(Vs,Ps,Es,varargin)
% Form a new random state around a uniform value (Given by M_InitUnfSt)
% VsOut = M_InitRndSt(Vs,Ps,Es,varargin)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'StNoise'))	
	Es.StNoise = Es.StSmall;
end;
Es.StNoise = [Es.StNoise(:)' 0 0];
 
% Prepare Uniform states if necessary
VsTmp = M_InitUnfSt(Vs,Ps,Es);

% Initilize randomization seed if relevant
if(Es.StNoise(2))
	rng(Es.StNoise(2));
end;

% Add a perturbation around the uniform state
VsOut = VsTmp + (rand(size(VsTmp))*2-1)*Es.StNoise(1); 

% If we know variables are positive, make sure they remain so	
if((isfield(Es,'NonNeg')) & (Es.NonNeg))
	VsOut = max(0,VsOut);
end;

end
