function VsNew=ChangeRes(Vs,Ps,Es,PsNew,varargin)
% Change a state in a current resolution to a different one
% VsNew=ChangeRes(Vs,Ps,Es,PsNew)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Is this a 1D system?
if((Ps.Nx==1) || (Ps.Ny==1))
    % Set correct resolution of the old system
	if(Ps.Nx==1)
		res=Ps.Ly/Ps.Ny;
	else	
		res=Ps.Lx/Ps.Nx;
	end;
    
    % Set correct resolution of the new system
	if(PsNew.Nx==1)
		resNew=PsNew.Ly/PsNew.Ny;
	else	
		resNew=PsNew.Lx/PsNew.Nx;
	end;
    
    VsNew = interp1(res*(0:Ps.Nx*Ps.Ny-1)',Vs,resNew*(1:PsNew.Nx*PsNew.Ny)'); 
    for ii=1:size(VsNew,2) 
        VsNew(isnan(VsNew(:,ii)),ii)=mean(Vs([1 end],ii));
    end;
    %VsNew = interp1(res*(1:Ps.Nx*Ps.Ny)',Vs,resNew*(1:PsNew.Nx*PsNew.Ny)'); 
    
else
    [xo,yo]=meshgrid((0:Ps.Ny-1)/(Ps.Ny-1)*Ps.Ly,(0:Ps.Nx-1)/(Ps.Nx-1)*Ps.Lx);
    [xi,yi]=meshgrid((1:PsNew.Ny)/PsNew.Ny*PsNew.Ly,(1:PsNew.Nx)/PsNew.Nx*PsNew.Lx);
    %size(yo),size(yi),size(reshape(Vs(:,1),Ps.Nx,Ps.Ny)')
    for ii=1:PsNew.VarNum
        ZI(:,:,ii) = interp2(xo,yo,reshape(Vs(:,ii),Ps.Nx,Ps.Ny),xi,yi);
    end;
    %size(ZI)
    VsNew = reshape(ZI,PsNew.Nx*PsNew.Ny,PsNew.VarNum);
     %warning('2D not supported');
end;

%VsNew(isnan(VsNew))=0;  % should be changed 

end

