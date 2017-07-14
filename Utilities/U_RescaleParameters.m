function [Vs,Ps,Es]=U_RescaleParameters(Vs,Ps,Es,varargin)
% A Utility function to recale parameters by a power of Es.RescaleFactor 
% [Vs,Ps,Es]=U_RescaleParameters(Vs,Ps,Es)
% Es.RescalePrm should have a list of parameters
% Es.RescalePow should have a list of powers to rescale by (default it 1)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% give default values
Es=InsertDefaultValues(Es,'RescaleFactor',[],'RescalePrm',[],'RescalePow',1);
% make sure these 2 fields are not empty
if(isempty(Es.RescalePrm) || isempty(Es.RescalePow)|| isempty(Es.RescaleFactor))
    error('Es.RescaleFactor, Es.RescalePrm and Es.RescalePow all need to be defined for rescaling.');
end; 

if(length(Es.RescalePow)<length(Es.RescalePrm))
    Es.RescalePow = [Es.RescalePow(:) ; repmat(Es.RescalePow(1),length(Es.RescalePrm),1)];
end;

oldnx = Ps.Nx;
oldny = Ps.Ny;

% load original values of parameters that will be rescaled
basevals=LoadParmList(Vs,Ps,Es,Es.RescalePrm);

for ii=1:length(Es.RescalePrm) % rescale each parameter
    tmpvals(ii)=basevals(ii)*(Es.RescaleFactor^Es.RescalePow(ii));
end;
% save these parameters into Ps/Es
[Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,tmpvals,Es.RescalePrm);

if(~(oldnx==Ps.Nx) || ~(oldny==Ps.Ny)) % change Vs if necessary
    [Vs,Ps,Es]=InitilizeState(mean(Vs,1),Ps,Es);
end;

end