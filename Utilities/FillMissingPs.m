function [Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es,varargin)
% A Utility function to fill in missing parameters
% [Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es)
% This function sets up an integration function if one is not defined,
% corrects the grid&system sizes for ill-defined 1D systems

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Make sure we have our 2 basic function handles (local, spatial)
if(~(isfield(Ps,'LocFunc')))
    error('The local model function Ps.LocFunc must be defined.');
end;
if(~(isfield(Ps,'SpaFunc')))
    error('The spatial model function Ps.SpaFunc must be defined.');
end;
% Make sure we have our boundary conditions and integration function set-up
if(~(isfield(Ps,'Bc')))
    Ps.Bc=0;
end;
if(~(isfield(Ps,'IntegFunc'))) 
    if(strcmp(func2str(Ps.SpaFunc),'S_RD')) % For reaction-diffusion
        if(Ps.Bc==0)
            Ps.IntegFunc=@I_PSRD; % Setup psedu-spectral if using periodic boundary conditions
        else
            Ps.IntegFunc=@I_FDSIMP; % otherwise use implicit method
        end;
    elseif(strcmp(func2str(Ps.SpaFunc),'S_LD'))
            Ps.IntegFunc=@I_FDSIMP; % For linear-derivatives, use implicit method
    else
        Ps.IntegFunc=@I_FDE;  % Assuming non linear-derivatives, use simple Euler method
    end;
end;
    
% Now look for grid&system size - check for 1D system not well defined
if(~(isfield(Ps,'Lx'))||~(isfield(Ps,'Nx'))) 
    if(~(isfield(Ps,'Ly'))||~(isfield(Ps,'Ny'))) 
        error('Either Lx&Nx or Ly&Ny must be defined in Ps.');
    else
        Ps.Lx=Ps.Ly; 
        Ps.Nx=Ps.Ny; 
        Ps.Ly=0;
        Ps.Ny=1;
    end;
elseif(~(isfield(Ps,'Ly'))||~(isfield(Ps,'Ny'))) 
	Ps.Ly=0;
	Ps.Ny=1;
end;
if(Ps.Nx<2) % Preferring system on the X axis, and not Y
    Ps.Lx=Ps.Ly; 
	Ps.Nx=Ps.Ny; 
	Ps.Ly=0;
	Ps.Ny=1;
end;
if(Ps.Ny<1)
    Ps.Ly=0;
	Ps.Ny=1;
end;

if(~(isfield(Ps,'VarNum')))
    if(strcmp(func2str(Ps.SpaFunc),'S_RD')) % For reaction-diffusion
        Ps.VarNum=length(Ps.Ds);
    else
        error('Variable number (Ps.VarNum) is undefined!');
    end;
end;

% TODO: Add sorting out the field order, if relevant

end


