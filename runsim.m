function [VsOut,dummyhist]=runsim(Vs,Ps,Es,varargin)
% Run integrator (Ps.IntegFunc) for a predefined time (Es.Tdest)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

Vs = Vs(:,:,1); % We're ignoring other states (other than the first one)

if(~isfield(Es,'Tmode'))
	Es.Tmode = 'none';
end;

% Calculate any matrices and other auxiliary data before run
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);

% Calculate time step automatically if relevant
if(strcmp(Es.Tmode,'auto'))
	Es.Tstep = EvaluateTS(Vs,Ps,Es);
end;

% Run integration
VsOut = Ps.IntegFunc(Vs,Ps,Es);		

dummyhist=[];
end