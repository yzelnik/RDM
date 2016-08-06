function goodts=EvaluateTS(Vs,Ps,Es,varargin)
% Evaluate a good time step to use with an integrator
% goodts=EvaluateTS(Vs,Ps,Es,IntFunc)
% Returns the timestep to use if successful, 0 otherwise
% IntFunc is the integrator function to use

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
Es.Tmode = 'none'; % Making sure this function does not call itself

stepnum = 100;  % how many steps to run each time
factor  = 2;    % muliply ts by how much each time?
mints   = 1e-7; % minimum value of ts to start from
maxts   = 1e+1; % max value of ts allowed
thresh  = 1e-1; % threshold of noise that is considered bad integration

% Start things by calculating integrating with slow time-step, and score=0
Es.Tstep=mints;
gs=runsim(Vs,Ps,Es,'Es.SkipWarning',1,'Es.Tdest',stepnum*Es.Tstep);
nrm = T_L2Norm(gs,Ps,Es);
score=0;

while (score<thresh) && (Es.Tstep<maxts)
	Es.Tstep=Es.Tstep*factor; % Increase time-step
    twofrms = runframes(Vs,Ps,Es,'Es.SkipWarning',1,'Es.Frames',[1/factor 1]*stepnum*Es.Tstep);
    score = T_L2Norm(gs - twofrms(:,:,1),Ps,Es);  % compare old run to new one
    %disp([Es.Tstep score])
    gs = twofrms(:,:,2);    % new run now becomes old

end;
%score
goodts = Es.Tstep/(factor^2);   % we found a "bad" time-step, no decrease back to "safe levels"

%tmp=IntFunc(Vs,Ps,Es,'Es.SkipWarning',1,'Es.Tstep',goodts,'Es.Tdest',stepnum*10*goodts);
%T_L2Norm(tmp,Ps,Es)
end



% OLD CODE - if becomes relevant

%firstTs=1;
%thresh=1e-9;
%smallTS=1e-20;
%Es.Tstep=firstTs;

% while flag==0
%Es.Tdest=runnum*Es.Tstep;
	%VsOut=IntFunc(Vs,Ps,Es,'Es.SkipWarning',1);
	
    %tot=sum(sum(abs(VsOut-Vs)));
    %[log10(tot) log10(Es.Tstep)]
	%if (isfinite(tot) && (tot<1/thresh) )
	%	flag=1;
	%end;
	%if(Es.Tstep<=smallTS)
	%	flag=2;
	%end;
%end;


%if flag==1
%	goodts=Es.Tstep;
%else
%	goodts=0;
%end;