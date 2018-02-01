function [Vs,Ps,Es]=U_RescaleParameters(Vs,Ps,Es,varargin)
% A Utility function to recale parameters by a power of Es.RescaleFactor 
% [Vs,Ps,Es]=U_RescaleParameters(Vs,Ps,Es)
% Es.RescalePrm should have a list of parameters
% Es.RescalePow should have a list of powers to rescale by (default it 1)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% give default values
Es=InsertDefaultValues(Es,'RescaleFactor',[],'RescalePrm',[],'RescalePow',1,'PrmInCell',[]);
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

tiny = 1e-10;
% try and round off rescaling value if relevant
if(mod(Es.RescaleFactor,1)<tiny || mod(-Es.RescaleFactor,1)<tiny)
    Es.RescaleFactor=round(Es.RescaleFactor);
elseif(mod(10*Es.RescaleFactor,1)<tiny || mod(-10*Es.RescaleFactor,1)<tiny)
    Es.RescaleFactor = round(10*Es.RescaleFactor)/10;
elseif(mod(100*Es.RescaleFactor,1)<tiny || mod(-100*Es.RescaleFactor,1)<tiny)
    Es.RescaleFactor = round(100*Es.RescaleFactor)/100;
end;
    

for ii=1:length(Es.RescalePrm) % rescale each parameter
    tmpvals(ii)=basevals(ii)*(Es.RescaleFactor^Es.RescalePow(ii));
end;

% save these parameters into Ps/Es
BackupPrmInCell=Es.PrmInCell; % In case we are using Es.PrmInCell, make sure we don't use it here
Es.PrmInCell=[];
[Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,tmpvals,Es.RescalePrm);
Es.PrmInCell=BackupPrmInCell; % In case we are using Es.PrmInCell, make sure we don't use it here

if(~(oldnx==Ps.Nx) || ~(oldny==Ps.Ny)) % change Vs if necessary
    Es.InitActive=0; % make sure we actually (re)initilize the state
    [Vs,Ps,Es]=InitilizeState(mean(Vs,1),Ps,Es);
end;

end

