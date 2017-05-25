function [Vs,Ps,Es]=U_SetupRndTime(Vs,Ps,Es,varargin)
% A Utility function to setup temporal heterogeneity in parameters 
% [Vs,Ps,Es]=U_SetupRndTime(Vs,Ps,Es)

error('function is currently not functional');

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% give default values
Es=InsertDefaultValues(Es,'RndTimePrm',[],'RndTimeVal',[]);
% make sure these 2 fields are not empty
if(isempty(Es.RndTimePrm) || isempty(Es.RndTimeVal))
    error('Both Es.RndTimePrm and Es.RndTimeVal need to be defined for setup.');
end; 

if(~iscell(Es.RndTimePrm)) % wrap in cell array if needed
    Es.RndTimePrm={Es.RndTimePrm};
end;

if(size(Es.RndTimeVal,1)<length(Es.RndTimePrm))
    error('missing info in Es.RndTimePrm');
end;

Es.RndTimeVal(1,size(Es.RndTimeVal,2)+3)=0; % adding zeros for buffer

% load average values of parameters that will be redefined
basevals=LoadParmList(Vs,Ps,Es,Es.RndTimePrm);

% go oer list of parameters, and calculate their random dist.
for ii=1:length(Es.RndTimePrm)
    if(Es.RndTimeVal(ii,2)==0) %uniform dist.
        %tmpvals{ii}=basevals(ii)+(rand(Ps.Nx*Ps.Ny,1)-0.5)*Es.RndTimeVal(ii,1);
    elseif(Es.RndTimeVal(ii,2)==1) % normal dist.
        %tmpvals{ii}=basevals(ii)+randn(Ps.Nx*Ps.Ny,1)*Es.RndTimeVal(ii,1);
    else
        error('undefined type of random distribution in Es.RndTimeVal');
    end;
end;
% save these parameters into Ps/Es
%[Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,tmpvals,Es.RndTimePrm);


end
