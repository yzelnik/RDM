function [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es,varargin)
% A Utility function to setup the spatial data (spatial matrix, flags, etc)
% [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es.SetupMode = 1; % Indicate we're in setup mode
Es.JacMode   = 0; % By default, we do not want a Jacobian
matlen = Ps.Nx*Ps.Ny*Ps.VarNum; % size of expected spatial matrix (if relevant)

% Try and setup the local function, if relevant
tmp = Ps.LocFunc(Vs,Ps,Es); 
% Update Ps, if indeed we got a struct back
if(isstruct(tmp)) Ps=tmp; end; 

% Erase previous flags if they exist
%Es.SpaMatUse = [];

% Try and setup the spatial matrix, and/or other such info of SpaFunc
tmp = Ps.SpaFunc(Vs,Ps,Es); 
% Update the Ps, if indeed we got a struct back
if(isstruct(tmp)) 
    Ps=tmp; 
else
    % If the size is right, assume we got a Spatial-Matrix back, so store it
    if(size(tmp,1)==matlen) && (size(tmp,2)==matlen)
        Ps.SpaMat=tmp;
    end;
end;

% Check if a spatial matrix exists and is not empty
if((isfield(Ps,'SpaMat')) && ~isempty(Ps.SpaMat))
	Es.SpaMatUse = 1;
else
    Es.SpaMatUse = 0;
end;

% Try and setup the integration function (without running it in time)
tmp = Wrap4Update(Vs,Ps,Es,Ps.IntegFunc,'Es.TimeDst',0);

% Update the Ps, if indeed we got a struct back
if(isstruct(tmp)) Ps=tmp; end; 

if((~isfield(Es,'JacNum')) || (isempty(Es.JacNum))) 
    locjac = Wrap4Update(Vs,Ps,Es,Ps.LocFunc,'Es.JacMode',1);
    if((size(locjac,1)==matlen) && (size(locjac,2)==matlen))
        Es.JacNum=0;
    else
        Es.JacNum=1;
    end;
end;
Es.SetupMode = 0;	% Back to normal

end
