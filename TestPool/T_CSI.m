function [flag,score]=T_CSI(Vs,Ps,Es,varargin)
% Check Stability by Integration of a numerical solution
% [flag,score]=T_CSI(Vs,Ps,Es)
% Returns 1 if the state is stable, 0 otherwise

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Initialization
Es.TimeMax = Es.TimeDst;

if(isfield(Es,'CsiThresh'))
    %disp(999);
    Es.SsThresh=Es.CsiThresh;
end;
% Integrate in time
VsOut=run2ss(Vs+rand(size(Vs))*(Es.StSmall/100),Ps,Es,'Es.OlDraw',0,'Es.TestFunc',[],'Es.NoWarning',1);

Es.SsThresh=Es.StSmall*10;

% Compare the new and old states
[flag,score]=CheckSS(cat(3,Vs,VsOut),Ps,Es);	
disp([score Es.SsThresh]);

end