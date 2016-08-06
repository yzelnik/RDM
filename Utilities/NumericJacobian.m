function [NumJac,AnaJac]=NumericJacobian(Vs,Ps,Es,varargin)
% Find an approxmation for a Jacobian numerically
% NumJac=NumericJacobian(Vs,Ps,Es,varargin)
% To compare the two jacobians, use as: [NumJac,AnaJac]=NumericJacobian(Vs,Ps,Es,varargin)

% Setup the spatial matrix and auxiliary flags
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);

% Initilize some basic values
szfactor = 1e2;
delta = max(eps*szfactor,Es.STsmall/szfactor);	% Make sure delta is a very small value, but not too small
syslen = Ps.Nx*Ps.Ny;
invdel = 1/delta;
Es.fmod=0;	% Just making sure we're not getting the analytical jacobian

% Now start things off. First calcualte the current right-hand-side of the pde
% After that, approximate derivatives, changing one point/site at a time

if(Es.UseSM)   % Get the right-hand-side of the PDE (no time derivatives)
    rhs0 = reshape(Ps.LocFunc(Vs,Ps,Es),syslen*Ps.Vnum,1) + Ps.SpaMat*Vs(:); 
else        % if we don't use SM (spatial matrix) than use the spatial function directly
    rhs0 = reshape( Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es) ,syslen*Ps.Vnum,1);
end; 

for ii=1:syslen	% Go over each data point
    for jj=1:Ps.Vnum
    	Vs(ii,jj)=Vs(ii,jj)+delta;	% Change the value at one point by a small amount (delta)
    	if Es.updateSM	
            Ps.SpaMat = Ps.SpaFunc(Vs,Ps,Es);  % Use this if the spatial matrix needs to be updated online
    	end;
    	if(Es.UseSM)   % Get the right-hand-side of the PDE (no time derivatives)
    	    rhs1 = reshape(Ps.LocFunc(Vs,Ps,Es),syslen*Ps.Vnum,1) + Ps.SpaMat*Vs(:); 
    	else        % if we don't use SM (spatial matrix) than use the spatial function directly
            rhs1 = reshape( Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es) ,syslen*Ps.Vnum,1);
    	end;  
    	NumJac(:,(jj-1)*syslen+ii) = (rhs1-rhs0)*invdel; % Calcualte one column of the jacobian matrix
    	Vs(ii,jj)=Vs(ii,jj)-delta;	% Change back the value at one point to its original one
    end;
end;

if(nargout>1) % Comparison of numerical jacobian to analytical one
    Es.Jcob = 1;	% Delete this line soon, old version (still needs change across the board)
    Es.fmod = 1;	% Request a jacobian

    AnaJac = Ps.LocFunc(Vs,Ps,Es)+Ps.SpaFunc(Vs,Ps,Es);
    
end;

end
