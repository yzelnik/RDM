function VsOut=M_InitSpots(Vs,Ps,Es,varargin)
% Mix two states to form a uniform state with a few spots
% VsOut=M_InitSpots(Vs,Ps,Es)
% Spots size and number according to: Es.InitPrm = [radius,spotnum]
MinSysLen = 10;

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'InitPrm') || isempty(Es.InitPrm))
	error('Need at least one value in Es.InitPrm (radius).');
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
Es.InitPrm
if(length(Es.InitPrm)<2)
    spotnum=1;
else
    spotnum=Es.InitPrm(2);
end;

% Initilize grid and mask
[yy,xx] = meshgrid(1:Ps.Ny,1:Ps.Nx);
mask = zeros(Ps.Nx,Ps.Ny);
temp = mask;

for ii=1:spotnum % Go over spots
    rx = rand(1)*Ps.Nx;
    ry = rand(1)*Ps.Ny;
    mask = (sqrt((xx-rx).^2+(yy-ry).^2)<=Es.InitPrm(1));
  
    temp(logical(mask(:)))=1;
end;

temp=repmat(reshape(temp,Ps.Nx*Ps.Ny,1),1,Ps.VarNum);

% Apply mask
VsOut = temp.*Vs(:,:,1) + (1-temp).*Vs(:,:,2);

end
