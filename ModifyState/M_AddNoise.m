function VsOut = M_AddNoise(Vs,Ps,Es,varargin)
% Add random noise to a given state using Es.StNoise (or Es.StSmall)
% VsOut = M_AddNoise(Vs,Ps,Es)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'StNoise'))	
	Es.StNoise = Es.StSmall;
end;
Es.StNoise = [Es.StNoise(:)' 0 0];
 
% Initilize randomization seed if relevant
if(Es.StNoise(2))
	rng(Es.StNoise(2));
end;

% Add a perturbation around the uniform state
VsOut = Vs + (rand(size(Vs))*2-1)*Es.StNoise(1); 

% If we know variables are positive, make sure they remain so	
if((isfield(Es,'NonNeg')) & (Es.NonNeg))
	VsOut = max(0,VsOut);
end;

end
