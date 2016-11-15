function testres=gettest(Vs,Ps,Es,varargin)
% Get test(s) results from Es.TestFunc on state(s) Vs
% st=gettest(Vs,Ps,Es)

% Default first extra input is for test function(s) to use
if(~mod(nargin,2)) varargin = ['Es.TestFunc' varargin]; end;

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Initilize state if necessary
[Vs,Ps,Es]=InitilizeState(Vs,Ps,Es);
% Setup the spatial matrix and auxiliary flags
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);	

% Find the steady-state
testres=PerformTests(Vs,Ps,Es,Es.TestFunc);

end

