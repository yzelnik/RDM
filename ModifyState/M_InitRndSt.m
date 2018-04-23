
function VsOut = M_InitRndSt(Vs,Ps,Es,varargin)
% Form a new random state around a uniform value (Given by M_InitUnfSt)
% VsOut = M_InitRndSt(Vs,Ps,Es,varargin)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'NonNeg',0,'StNoise',Es.StSmall);
Es.StNoise = [Es.StNoise(:)' 0 0];
 
% Prepare Uniform states if necessary
VsTmp = M_InitUnfSt(Vs,Ps,Es);

% Initilize randomization seed if relevant
if(Es.StNoise(2))
	rng(Es.StNoise(2));
end;

% Add a perturbation around the uniform state
VsOut = VsTmp + (rand(size(VsTmp))*2-1)*Es.StNoise(1); 

if(Es.NonNeg) % make sure values are not negative, if relevant
	VsOut = max(0,VsOut);
end;

end
