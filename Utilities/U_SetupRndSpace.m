function [Vs,Ps,Es]=U_SetupRndSpace(Vs,Ps,Es,varargin)
% A Utility function to setup spatial heterogeneity in parameters 
% [Vs,Ps,Es]=U_SetupRndSpace(Vs,Ps,Es)
% Es.RndSpacePrm should have a list of parameters
% Es.RndSpaceVal should have a list of values
% where the first column shows the width of distribution
% and (a possible) second column the type of distribution
% (0/none = uniform, 1=normal)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% give default values
Es=InsertDefaultValues(Es,'RndSpacePrm',[],'RndSpaceVal',[]);
% make sure these 2 fields are not empty
if(isempty(Es.RndSpacePrm) || isempty(Es.RndSpaceVal))
    error('Both Es.RndSpacePrm and Es.RndSpaceVal need to be defined for setup.');
end; 

if(~iscell(Es.RndSpacePrm)) % wrap in cell array if needed
    Es.RndSpacePrm={Es.RndSpacePrm};
end;

if(size(Es.RndSpaceVal,1)<length(Es.RndSpacePrm))
    error('missing info in Es.RndSpacePrm');
end;

Es.RndSpaceVal(1,size(Es.RndSpaceVal,2)+3)=0; % adding zeros for buffer

% load average values of parameters that will be redefined
basevals=LoadParmList(Vs,Ps,Es,Es.RndSpacePrm);

% go oer list of parameters, and calculate their random dist.
for ii=1:length(Es.RndSpacePrm)
    if(Es.RndSpaceVal(ii,2)==0) %uniform dist.
        tmpvals{ii}=basevals(ii)+(rand(Ps.Nx*Ps.Ny,1)-0.5)*Es.RndSpaceVal(ii,1);
    elseif(Es.RndSpaceVal(ii,2)==1) % normal dist.
        tmpvals{ii}=basevals(ii)+randn(Ps.Nx*Ps.Ny,1)*Es.RndSpaceVal(ii,1);
    else
        error('undefined type of random distribution in Es.RndSpaceVal');
    end;
end;
% save these parameters into Ps/Es
[Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,tmpvals,Es.RndSpacePrm);


end
