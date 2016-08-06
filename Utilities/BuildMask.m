function mask=BuildMask(parms,Ps,varargin)
% Build a mask in size specified by Ps, according to parms
% mask=BuildMask(parms,Ps)
% Paramaters can be given in relative (0..1) terms, or in Nx/Ny terms
% If any value is negative, an inverse mask (1-mask) is used
% The parameters in parms are interpeted by the number of them specified:
% 1 - the system is cut into two parts, in this location (in X)
% 2 - A rectangle from (0,0) to (X,Y), with parms = [X,Y,R]
% 3 - A circle with center of (X,Y) and radius R, with parms = [X,Y,R]
% 4 - A rectangle from (X1,Y1) to (X2,Y2), with parms = [X1,Y1,X2,Y2]

if(~isstruct(Ps))	% Build Ps struct online, mostly as a shortcut
	temp = [Ps(:); 1];
	Ps = struct('Nx',temp(1),'Ny',temp(2));
end;

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters([],Ps,[],varargin{:});

if(~length(parms))
	parms = [0.5];
end;

if(sum(parms<0)>0)
	parms = abs(parms);
	negflag = 1;
else
	negflag = 0;
end;

% First check if this is a 1D system
if((Ps.Nx==1) || (Ps.Ny==1))
	Ps.Nx=Ps.Nx*Ps.Ny;
	if(sum(abs(parms)<=1)>0)	% If paramteres are given in relative terms, give them in Nx dimensions
		parms = round(parms(:).*Ps.Nx);
	end;
	
	if(length(parms)==1)
		mask = [ones(parms(1),1) ; zeros(Ps.Nx-parms(1),1)];
	else
		mask = zeros(Ps.Nx,1);
		mask(max(1,parms(1)-parms(2)):min(Ps.Nx,parms(1)+parms(2)))=1;
	end;

else  % Assuming this is a 2D system
	if(sum(abs(parms)<=1)>0)	% If paramteres are given in relative terms, give them in Nx&Ny dimensions
		temp = repmat([Ps.Nx ;Ps.Ny],4,1);
		parms = round(parms(:).*temp(1:length(parms)));
	end;
	
	if(length(parms)==1)
		mask = [ones(parms(1),Ps.Ny) ; zeros(Ps.Nx-parms(1),Ps.Ny)]';
		
	elseif (length(parms)==2)
		mask = zeros(Ps.Ny,Ps.Nx);
		mask(1:parms(2),1:parms(1))=1;
	elseif(length(parms)==3)
		mask = zeros(Ps.Nx,Ps.Ny);
		[xx yy] = meshgrid(1:Ps.Nx,1:Ps.Ny);
		mask = (sqrt((xx-parms(1)).^2+(yy-parms(2)).^2)<=parms(3));
	else
		mask = zeros(Ps.Ny,Ps.Nx);
		mask(parms(2):parms(4),parms(1):parms(3))=1;
	end;
	
	mask = reshape(mask',length(mask(:)),1);

end

if(negflag)
	mask = 1-mask;
end;

end
