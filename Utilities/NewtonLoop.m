function [VsOut,ind,residual]=NewtonLoop(Vs,Ps,Es,varargin)
% Repeat a Newton loop until threshold is reached
% [VsOut,ind,residual]=NewtonLoop(Vs,Ps,Es)
% Ps.LocFunc & Ps.SpaFunc are the functions for the local & non-local parts
% The threshold is set by Es.SSthresh (or 1e-10 by default)
maxloop = 20; % max number of iterations

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Setup the spatial matrix and auxiliary flags
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es,varargin{:});


%if(~isfield(Es,'fmod'))
%   Es.fmod=0;
%end;
posflag = 0;
if((isfield(Es,'posflag')) & (Es.posflag))
    posflag = 1;
end;
if(~isfield(Es,'SSthresh'))
   Es.SSthresh=1e-10;
end;
if(~isfield(Es,'NumJac') || Es.NumJac==0)
    NumJac=0;
else
    NumJac=1;
end;
if(~isfield(Es,'NewtonRate') || Es.NewtonRate==0)   % For slower convegence
    rate = 1;
else
    rate = Es.NewtonRate;
end;
syslen = Ps.Nx * Ps.Ny;
totlen = syslen * Ps.Vnum;


ind = 1;
residual = inf;
while (ind<maxloop) && (residual>Es.SSthresh)
    if(Es.UseSM)   % Get the right-hand-side of the PDE (no time derivatives)
        rhs = Ps.LocFunc(Vs,Ps,Es) + reshape(Ps.SpaMat*Vs(:),syslen,Ps.Vnum); 
    else        % if we don't use SM (spatial matrix) than use the spatial function directly
        rhs = Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es);
    end; 
    residual = (sqrt(sum(rhs(:).^2))/totlen);
    
    if(NumJac==1)
        jac = NumericJacobian(Vs,Ps,Es);
        %disp('numeric!');
    else
        Es.fmod = 1;	% Request a jacobian
        jac = Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es);
        Es.fmod = 0;	% Back to normal
    end;
    % Use the inverse of the jac matrix on the rhs vector to update the state vector
    VsChange = jac \ reshape(rhs,totlen,1); 
    Vs = Vs - rate*reshape(VsChange,syslen,Ps.Vnum);

    if Es.updateSM
        Ps.SpaMat = Ps.SpaFunc(Vs,Ps,Es);  % Use this if the spatial matrix needs to be updated online
    end;
    if posflag  % consider deleting
        	Vs = max(0,Vs);
    end;
    ind = ind+1;
    %history(:,ind)=Vs; % delete this line
end;
%plot(history)% delete this line
VsOut = Vs;
%disp([ind log10(residual)])% delete this line

end
