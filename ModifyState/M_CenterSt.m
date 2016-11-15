function VsOut=M_CenterSt(Vs,Ps,Es,varargin)
% Center a given state for variable Es.VarInd(1) (=1 by default)
% This works by centerning on the maximal region of the variable

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% VarInd -> Main variable to work on
if isfield(Es,'VarInd')                  
    Es.VarInd = Es.VarInd(1);
else
    Es.VarInd = 1;
end;

% state to work on
st = Vs(:,Es.VarInd,1);

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
    VsOut = M_ShiftSt(Vs,Ps,Es,'Es.ShiftPrm',round(-mxloc));
 
else  % Assuming this is a 2D system
    % create two upside-down-v-shapes, for x&y
    maxcenterx = reshape(repmat([0.5:Ps.Nx/2 Ps.Nx/2:-1:0.5]',1,Ps.Ny),Ps.Nx*Ps.Ny,1)';
    maxcentery = reshape(repmat([0.5:Ps.Nx/2 Ps.Nx/2:-1:0.5],Ps.Nx,1),Ps.Nx*Ps.Ny,1)';
    
    % get score for each cyclic shift, for both x&y
    for ii=1:Ps.Nx 
        tmp = M_ShiftSt(Vs,Ps,Es,'Es.ShiftPrm',[ii 0]);
        midx(ii)=maxcenterx*tmp(:,Es.VarInd);
    end; 
    for ii=1:Ps.Nx 
        tmp = M_ShiftSt(Vs,Ps,Es,'Es.ShiftPrm',[0 ii]);
        midy(ii)=maxcentery*tmp(:,Es.VarInd);
    end; 
    % find the location of the best score
    [~,mxlocx]=max(midx);
    [~,mxlocy]=max(midy);
    
    % center on mxlocx,mxlocy
    VsOut = M_ShiftSt(Vs,Ps,Es,'Es.ShiftPrm',[mxlocx mxlocy]);
end;


end
