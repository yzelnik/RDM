function VsOut = InitStateByMask(Vs,Ps,Es,varargin)
% Initiate a state using initial values given in Vs, and a mask given by Es.StMask
% VsOut = InitStateByMask(Vs,Ps,Es,varargin)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

Vs = UpdateInitUniform(Vs,Ps,Es);
mask = Es.StMask(:);

% Build a matrix (of size pnum*vnum) with 1's where the relevant initial values should be used
base = zeros(length(mask),Ps.Vnum);
base(sub2ind(size(base),(1:length(mask))',mask)) = 1;

% Build the inital state
VsOut = base * Vs;


end
