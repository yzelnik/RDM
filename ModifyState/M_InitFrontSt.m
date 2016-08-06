function VsOut=M_InitFrontSt(Vs,Ps,Es,varargin)
% Mix two states to form a front
% Es.InitParms gives the relative part(s) of the two domains 
% VsOut=M_InitFrontSt(Vs,Ps,Es)
MinSysLen = 10;

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'InitParms'))
	Es.InitParms = [0.5];
end;

% Make mask per variable
mask=[];
for ind=1:length(Es.InitParms)
    tmpmask = BuildMask(Es.InitParms(ind),Ps);
    mask = [mask tmpmask];
end;
% If InitParms is not long enough, repeat the mask as necessary
if(size(mask,2)<Ps.Vnum)
    mask = [mask repmat(mask,1,Ps.Vnum-size(mask,2))];
end;

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
