function VsOut=M_CenterSt(Vs,Ps,Es,varargin)
% Center a given state for variable Es.Vind(1) (=1 by default)
% This works by centerning on the maximal region of the variable

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Vind -> Main variable to work on
if isfield(Es,'Vind')                  
    Es.Vind = Es.Vind(1);
else
    Es.Vind = 1;
end;

% state to work on
st = Vs(:,Es.Vind,1);

% First check if this is a 1D system
if((Ps.Nx==1) || (Ps.Ny==1))
    % 1D case
	Ps.Nx=Ps.Nx*Ps.Ny;
    % create an upside-down-v-shape
    maxcenter=[0.5:Ps.Nx/2 Ps.Nx/2:-1:0.5];
    
    % get score for each cyclic shift
    for ii=1:Ps.Nx 
        midx(ii)=maxcenter*circshift(st(:,1),ii); 
    end; 
    % find the best score
    [~,mxloc]=max(midx);
    % center on mxloc
    VsOut = M_ShiftSt(Vs,Ps,Es,'Es.ShiftParms',round(-mxloc));
 
else  % Assuming this is a 2D system
    % create two upside-down-v-shapes, for x&y
    maxcenterx = reshape(repmat([0.5:Ps.Nx/2 Ps.Nx/2:-1:0.5]',1,Ps.Ny),Ps.Nx*Ps.Ny,1)';
    maxcentery = reshape(repmat([0.5:Ps.Nx/2 Ps.Nx/2:-1:0.5],Ps.Nx,1),Ps.Nx*Ps.Ny,1)';
    
    % get score for each cyclic shift, for both x&y
    for ii=1:Ps.Nx 
        tmp = M_ShiftSt(Vs,Ps,Es,'Es.ShiftParms',[ii 0]);
        midx(ii)=maxcenterx*tmp(:,Es.Vind);
    end; 
    for ii=1:Ps.Nx 
        tmp = M_ShiftSt(Vs,Ps,Es,'Es.ShiftParms',[0 ii]);
        midy(ii)=maxcentery*tmp(:,Es.Vind);
    end; 
    % find the location of the best score
    [~,mxlocx]=max(midx);
    [~,mxlocy]=max(midy);
    
    % center on mxlocx,mxlocy
    VsOut = M_ShiftSt(Vs,Ps,Es,'Es.ShiftParms',[mxlocx mxlocy]);
end;


end
