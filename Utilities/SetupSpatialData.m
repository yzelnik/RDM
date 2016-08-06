function [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es,varargin)
% A Utility function to setup the spatial data (spatial matrix, flags, etc)
% [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Erase previous flags if they exist
Es.UseSM    = [];
Es.updateSM = [];

% Try and setup the spatial matrix, and/or other such info of SpaFunc
Es.SetupMode = 1;	
tmp = Ps.SpaFunc(Vs,Ps,Es); 
% Update the Ps, if indeed we got a struct back
if(isstruct(tmp)) Ps=tmp; end; 

% Check if a spatial matrix exists and it has the right size
if((isfield(Ps,'SpaMat')) && (size(Ps.SpaMat,1)==size(Ps.SpaMat,2)))
    Es.UseSM = 1;
else
    Es.UseSM = 0;
end;

% Try and setup initial stuff in the integration function
tmp = Ps.IntegFunc(Vs,Ps,Es,'Es.Tdest',0); 
% Update the Ps, if indeed we got a struct back
if(isstruct(tmp)) Ps=tmp; end; 

Es.SetupMode = 0;	% Back to normal


% Check if the Spatial matrix needs to be updated (due to Non-Linear Derivatives)
Es.updateSM = 0;	
if((isfield(Ps,'NLD')) & (~isempty(Ps.NLD)))
    Es.updateSM = 1;
end;



end
