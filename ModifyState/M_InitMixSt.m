function VsOut=M_InitMixSt(Vs,Ps,Es,varargin)
% Mix two states to form a new initial conditions state
% Es.InitPrm should contain a full mask for this, or the parameters to form the mask
% VsOut=M_InitMixSt(Vs,Ps,Es)
MinSysLen = 10;

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'InitPrm'))
	Es.InitPrm = [];
end;

% If InitPrm is long enough, assume its a full ready mask, otherwise build a mask
if(length(Es.InitPrm)>MinSysLen)
	mask = Es.InitPrm;	
else
	mask = BuildMask(Es.InitPrm,Ps);
end;
mask = repmat(mask,1,Ps.VarNum);

% If only mean-field values are given for Vs
if(size(Vs,1)<MinSysLen)
	Vs = M_InitUnfSt(Vs,Ps,Es);
end;

% If only one full base state is given, use a uniform as the second one
if(size(Vs,3)==1)
	Vs2 = M_InitUnfSt(mean(Vs,1),Ps,Es);
	Vs = cat(3,Vs,Vs2);
end;

% Apply mask
VsOut = mask.*Vs(:,:,1) + (1-mask).*Vs(:,:,2);

end
