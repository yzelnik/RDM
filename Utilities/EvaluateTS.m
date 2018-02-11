function goodts=EvaluateTS(Vs,Ps,Es,varargin)
% Evaluate a good time step to use with an integrator
% goodts=EvaluateTS(Vs,Ps,Es,IntFunc)
% Returns the timestep to use if successful, 0 otherwise
% IntFunc is the integrator function to use

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'Verbose',0,'TsMax',1e+1,'TsMin',1e-7);

Es.OlDraw = 0; % making sure we're not plotting anything...
Es.TsMode = 'none'; % Making sure this function does not call itself
% Initilize state if necessary
[Vs,Ps,Es]=InitilizeState(Vs,Ps,Es);

stepnum = 100;      % how many steps to run each time
factor  = 2;        % muliply ts by how much each time?
mints   = Es.TsMin; % minimum value of ts to start from
maxts   = Es.TsMax; % max value of ts allowed
thresh  = 1e-1;     % threshold of noise that is considered bad integration

% Start things by calculating integrating with slow time-step, and score=0
Es.TsSize=mints;
gs=runsim(Vs,Ps,Es,'Es.NoWarning',1,'Es.TimeDst',stepnum*Es.TsSize);

score=0;
while (score<thresh) && (Es.TsSize<maxts)
	Es.TsSize=Es.TsSize*factor; % Increase time-step
    twofrms = runframes(Vs,Ps,Es,'Es.NoWarning',1,'Es.Frames',[1/factor 1]*stepnum*Es.TsSize,'Es.FramesChoice',1:2,'Es.RecurFrames',[],'Es.DynPrm',[]);
    score = T_L2Norm(gs - twofrms(:,:,1),Ps,Es);  % compare old run to new one
    %disp([Es.TsSize score])
    gs = twofrms(:,:,2);    % new run now becomes old
end;
%score
goodts = Es.TsSize/(factor^2);   % we found a "bad" time-step, now decrease back to "safe levels"
% round it down into a nice number
tmp=10.^floor(log10(abs(goodts)));
goodts=floor(goodts/tmp)*tmp;

if(Es.Verbose)
    disp(sprintf('ts = %e',goodts));
end;

end



