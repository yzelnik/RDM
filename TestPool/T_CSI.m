function [flag,score]=T_CSI(Vs,Ps,Es,varargin)
% Check Stability by Integration of a numerical solution
% [flag,score]=T_CSI(Vs,Ps,Es)
% Returns 1 if the state is stable, 0 otherwise

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Initialization
Es.Tmax = Es.Tdest;

% Integrate in time
VsOut=run2ss(Vs+rand(size(Vs))*(Es.STsmall/100),Ps,Es,'Es.OLdraw',0);

Es.SSthresh=Es.STsmall*10;

% Compare the new and old states
[flag,score]=CheckSS(cat(3,Vs,VsOut),Ps,Es);	
disp([score Es.SSthresh]);

end