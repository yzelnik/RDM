function VsOut=MixStates(Vs,Ps,Es,varargin)
% Mix two or more states
% VsOut=MixStates(Vs,Ps,Es)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if(~isfield(Es,'MixParms'))
	Es.MixParms = [0.5];
end;

if(sum(Es.MixParms<0)>0)
	Es.MixParms = abs(Es.MixParms);
	NegMaskFlag = 1;
else
	NegMaskFlag = 0;
end;

% First check if this is a 1D system
if((Ps.Nx==1) || (Ps.Ny==1))
	Ps.Nx=Ps.Nx*Ps.Ny;
	if(sum(abs(Es.MixParms)<=1)>0)	% If paramteres are given in relative terms, give them in Nx dimensions
		Es.MixParms = round(Es.MixParms(:).*Ps.Nx);
	end;
	
	if(length(Es.MixParms)==1)
		mask = [ones(Es.MixParms(1),1) ; zeros(Ps.Nx-Es.MixParms(1),1)];
	else
		mask = zeros(Ps.Nx,1);
		mask(max(1,Es.MixParms(1)-Es.MixParms(2)):min(Ps.Nx,Es.MixParms(1)+Es.MixParms(2)))=1;
	end;

else  % Assuming this is a 2D system
	if(sum(abs(Es.MixParms)<=1)>0)	% If paramteres are given in relative terms, give them in Nx&Ny dimensions
		temp = repmat([Ps.Nx ;Ps.Ny],4,1);
		Es.MixParms = round(Es.MixParms(:).*temp(1:length(Es.MixParms)));
	end;
	
	if(length(Es.MixParms)==1)
		mask = [ones(Es.MixParms(1),Ps.Ny) ; zeros(Ps.Nx-Es.MixParms(1),Ps.Ny)]';
		
	elseif (length(Es.MixParms)==2)
		mask = zeros(Ps.Ny,Ps.Nx);
		mask(1:Es.MixParms(2),1:Es.MixParms(1))=1;
	elseif(length(Es.MixParms)==3)
		mask = zeros(Ps.Nx,Ps.Ny);
		[xx yy] = meshgrid(1:Ps.Nx,1:Ps.Ny);
		mask = (sqrt((xx-Es.MixParms(1)).^2+(yy-Es.MixParms(2)).^2)<=Es.MixParms(3));
	else
		mask = zeros(Ps.Ny,Ps.Nx);
		mask(Es.MixParms(2):Es.MixParms(4),Es.MixParms(1):Es.MixParms(3))=1;
	end;
	mask = repmat(reshape(mask',length(mask(:)),1),1,Ps.Vnum);

end

if(NegMaskFlag)
	mask = 1-mask;
end;

if(size(Vs,3)==1)
	bla = UpdateInitUniform(mean(Vs,1),Ps,Es);
	Vs2 = repmat(bla,Ps.Nx*Ps.Ny,1);
	Vs = cat(3,Vs,Vs2);
end;
	VsOut = mask.*Vs(:,:,1) + (1-mask).*Vs(:,:,2);
end
