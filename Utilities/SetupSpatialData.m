function [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es,varargin)
% A Utility function to setup the spatial data (spatial matrix, flags, etc)
% [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es.SetupMode = 1;

% Try and setup the local function, if relevant
tmp = Ps.LocFunc(Vs,Ps,Es); 
% Update the Ps, if indeed we got a struct back
if(isstruct(tmp)) Ps=tmp; end; 

% Erase previous flags if they exist
Es.SmUse    = [];
Es.SmUpdate = [];

% Try and setup the spatial matrix, and/or other such info of SpaFunc
tmp = Ps.SpaFunc(Vs,Ps,Es); 
% Update the Ps, if indeed we got a struct back
if(isstruct(tmp)) Ps=tmp; end; 

matlen = Ps.Nx*Ps.Ny*Ps.VarNum;

% Check if a spatial matrix exists and it has the right size
if((isfield(Ps,'SpaMat')) && (size(Ps.SpaMat,1)==matlen) && (size(Ps.SpaMat,2)==matlen))
    Es.SmUse = 1;
else
    Es.SmUse = 0;
end;

% Try and setup initial stuff in the integration function
EsTemp=Es;
EsTemp.TimeDst=0;
tmp = Ps.IntegFunc(Vs,Ps,EsTemp); 

% Update the Ps, if indeed we got a struct back
if(isstruct(tmp)) Ps=tmp; end; 

% Check if the Spatial matrix needs to be updated (due to Non-Linear Derivatives)
Es.SmUpdate = 0;	
if((isfield(Ps,'Nld')) && (~isempty(Ps.Nld)))
    Es.SmUpdate = 1;
end;

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
