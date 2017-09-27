
function VsOut=M_ShiftSt(Vs,Ps,Es,varargin)
% Shift a given state by Es.ShiftPrm
% Es.ShiftPrm can be given in relative (0..1) terms, or in Nx/Ny terms

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% First check if this is a 1D system
if((Ps.Nx==1) || (Ps.Ny==1))
	Ps.Nx=Ps.Nx*Ps.Ny;
    	shift=Es.ShiftPrm(1);
	if((abs(shift)<1)&&(abs(shift)>0))	% If paramtere is given in relative terms, give them in Nx dimensions
		shift = round(shift.*Ps.Nx);
	end;
    shift = mod(shift,Ps.Nx);
 
    VsOut = cat(1,Vs(shift+1:end,:,:),Vs(1:shift,:,:));
else  % Assuming this is a 2D system
    shiftx=Es.ShiftPrm(1);
    shifty=Es.ShiftPrm(2);
	if(((abs(shiftx)<1)&&(abs(shiftx)>0))||((abs(shifty)<1)&&(abs(shifty)>0)))	% If paramteres are given in relative terms, give them in Nx&Ny dimensions
        shiftx = round(shiftx*Ps.Nx);
        shifty = round(shifty*Ps.Ny);
	end;  
    shiftx = mod(shiftx,Ps.Nx);
    shifty = mod(shifty,Ps.Ny);
    
    temp = reshape(Vs,Ps.Nx,Ps.Ny,size(Vs,2),size(Vs,3));
    temp = cat(1,temp(shiftx+1:end,:,:,:),temp(1:shiftx,:,:,:));
    temp = cat(2,temp(:,shifty+1:end,:,:),temp(:,1:shifty,:,:));
    VsOut= reshape(temp,Ps.Nx*Ps.Ny,size(Vs,2),size(Vs,3));
    %VsOut = cat(1,Vs(Es.ShiftPrm(1):end,:,:),Vs(1:Es.ShiftPrm(1)-1,:,:));
    
end;


end
