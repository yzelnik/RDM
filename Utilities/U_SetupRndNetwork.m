function [Vs,Ps,Es]=U_SetupRndNetwork(Vs,Ps,Es,varargin)
% A Utility function to setup spatial heterogeneity via a network
% [Vs,Ps,Es]=U_SetupRndNetwork(Vs,Ps,Es)
% Es.RndNetFunc should have function to create a network
% Es.RndNetPrm should have a list of parameters for this function

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% give default values
Es=InsertDefaultValues(Es,'RndNetFunc',[],'RndNetPrm',[]);
% make sure these 2 fields are not empty
if(isempty(Es.RndNetFunc) || isempty(Es.RndNetPrm))
    error('Both Es.RndNetFunc and Es.RndNetPrm need to be defined for setup.');
end; 

Ps.Net = Es.RndNetFunc(Ps.Nx,Es.RndNetPrm);

end
