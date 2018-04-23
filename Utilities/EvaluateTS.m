function goodts=EvaluateTS(Vs,Ps,Es,varargin)
% Evaluate a good time step to use with an integrator
% goodts=EvaluateTS(Vs,Ps,Es)
% Uses  Es.TsMin and Es.TsMax for limits on possible time-step sizes (def=[1e-7 1e-1])
% And Es.EvalTsPrm = [stepnum factor thresh] for the iteration process (def=[100 2 0.1])
% where stepnum is the numer if integration steps each time, factor is how
% much to increase (=multiply) the time-step each iteration, and thresh defines bad integration 

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'TsMax',1e+1,'TsMin',1e-7,'EvalTsPrm',[100 2 0.1]);

Es.OlDraw = 0;      % Making sure we're not plotting anything
Es.TsMode = 'none'; % Making sure this function does not call itself
Es.RecurFrames =[]; % Simplyfing the runframes that we'll use
Es.DynPrm = [];     % Simplyfing the runframes that we'll use
Es.TestFunc = [];   % Simplyfing the runframes that we'll use

% Initilize state if necessary
[Vs,Ps,Es]=InitilizeState(Vs,Ps,Es);

mints   = Es.TsMin; % Minimum value of ts to start from (def=1e-7)
maxts   = Es.TsMax; % Max value of ts allowed (def=1e-1)
stepnum = Es.EvalTsPrm(1);	% How many steps to run each time (def=100)
factor  = Es.EvalTsPrm(2);	% Muliply ts by how much each time? (def=2)
thresh  = Es.EvalTsPrm(3);	% Threshold of noise that is considered bad integration (def=0.1)

% Start things by calculating integrating with slow time-step, and score=0
Es.TsSize=mints;
gs=runsim(Vs,Ps,Es,'Es.NoWarning',1,'Es.TimeDst',stepnum*Es.TsSize);

score=0;
while (score<thresh) && (Es.TsSize<maxts)
	Es.TsSize=Es.TsSize*factor; % Increase time-step
    twofrms = runframes(Vs,Ps,Es,'Es.NoWarning',1,'Es.Frames',[1/factor 1]*stepnum*Es.TsSize,'Es.FramesChoice',1:2);
    score = T_L2Norm(gs - twofrms(:,:,1),Ps,Es);  % Compare old run to new one

    gs = twofrms(:,:,2);    % New run now becomes old
end;

goodts = Es.TsSize/(factor^2);   % We found a "bad" time-step, now decrease back to "safe levels"

% Round it down into a nice number
tmp=10.^floor(log10(abs(goodts)));
goodts=floor(goodts/tmp)*tmp;

end



