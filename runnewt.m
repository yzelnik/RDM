function [VsOut,residual,ind]=runnewt(Vs,Ps,Es,varargin)
% Run Newton loop until threshold is reached (Es.SsThresh or 1e-10 by default)
% [VsOut,ind,residual]=runnewt(Vs,Ps,Es)

% Default first extra input is for the threshold to stop simulation
if(~mod(nargin,2)) varargin = ['Es.SsThresh' varargin]; end;

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Initilize state if necessary
[Vs,Ps,Es]=InitilizeState(Vs,Ps,Es);
% Calculate any matrices and other auxiliary data before run
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);

Vs = Vs(:,:,1); % We're ignoring other states (other than the first one)

% Run the Newton-loop
[VsOut,residual,ind]=NewtonLoop(Vs,Ps,Es);	

end
