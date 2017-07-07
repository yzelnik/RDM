function [Vs,Ps,Es]=U_SetupPoissonProcess(Vs,Ps,Es,varargin)
% A Utility function to setup A Poisson point process for runframes
% Es.PppPrm(1) is the average time between events
% [Vs,Ps,Es]=U_SetupPoissonProcess(Vs,Ps,Es)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% give default values
Es=InsertDefaultValues(Es,'PppPrm',[]);
% make sure these 2 fields are not empty
if(isempty(Es.PppPrm))
    error('Both Es.PppPrm needs to be defined for setup (including at least the rate of events).');
end; 

Es.PppPrm = [Es.PppPrm(:)' 0 0]; % padding zeros for buffer

res     = diff(Es.Frames(1:2));
avgjmp  = Es.PppPrm(1);
tottime = Es.Frames(end);

if(sum(abs(diff(Es.Frames)-res)>Es.TsSize)) % making sure Es.Frames is setup properly
    error('Es.Frames contains non-ordered frames');
end;

if(Es.PppPrm(2)==0) % regular Poisson process (no correlations)
    tmptms = cumsum(-avgjmp*log(rand(round(2*tottime/avgjmp),1)));
    Es.RecurFrames = unique([1; ceil(tmptms/res)]);
else
    error('correlated events not implemented yet.');
end;

end
