function [VsOut,dummyhist]=runsim(Vs,Ps,Es,varargin)
% Run integrator (Ps.IntegFunc) for a predefined time (Es.TimeDst)
% VsOut=runsim(Vs,Ps,Es)

% Default first extra input is for the time to run the simulation
if(~mod(nargin,2)) varargin = ['Es.TimeDst' varargin]; end;
    
% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'TsMode','none');
% Initilize state if necessary
[Vs,Ps,Es]=InitilizeState(Vs,Ps,Es);

% Calculate time step automatically if relevant
[Vs,Ps,Es]=SetupTimeStep(Vs,Ps,Es);
% Calculate any matrices and other auxiliary data before run
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);

Vs = Vs(:,:,1); % We're ignoring other states (other than the first one)

% Run integration
VsOut = Ps.IntegFunc(Vs,Ps,Es);		

dummyhist=[];
end