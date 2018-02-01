function [VsOut,history,time]=run2ss(Vs,Ps,Es,varargin)
% Run integrator until steady state is reached 
% [VsOut,history,time]=run2ss(Vs,Ps,Es)
% steady-state is set by ES.SsThresh (or, if not defined, by Es.StSmall) 

% Default first extra input is for the threshold to stop simulation
if(~mod(nargin,2)) varargin = ['Es.SsThresh' varargin]; end;
    
% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'TimeMax',inf,'TimeMin',0,'TsNum',1e2,'OlDraw',0,'TestFunc',[],'TsMode','none','NoWarning',0,'PlotFunc',@plotst);
% Initilize state if necessary
[Vs,Ps,Es]=InitilizeState(Vs,Ps,Es);

% Unless otherwise specified, check if in one "time unit" a threshold of StSmall is passed
if(~isfield(Es,'SsThresh'))	
	Es.SsThresh = Es.StSmall/10;
end;
% In case we want to run several test functions
if(iscell(Es.TestFunc))  
    Es.TestList=Es.TestFunc;
    Es.TestFunc=@T_MultiTest;
end;

% Calculate time step automatically if relevant
[Vs,Ps,Es]=SetupTimeStep(Vs,Ps,Es);
% Calculate any matrices and other auxiliary data before run
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);

Vs = Vs(:,:,1); % Take only first state, if there are more than one.
Es.TimeDst = Es.TsSize*Es.TsNum;

% Set up GUI for online drawing
FlagStop = 0;
if Es.OlDraw
    clf;
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.7 0.09,0.1],'string','Stop','callback',{@stopb});
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.5 0.09,0.1],'string','Pause','callback',{@pauseb});
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.3 0.09,0.1],'string','Finish','callback',{@finishb});

    if(isempty(Es.TestFunc)) || (Es.OlDraw==1)  % If no test was defined, but we are using Online-draw, make something up...
        testtext='log10 of difference';
    else
        testtext=regexprep(func2str(Es.TestFunc),'_','-');
    end;
    
end

% Recalculating Es.SsThresh for step size...
Es.SsThresh=Es.SsThresh*Es.TimeDst;
%disp([ Es.TimeDst Es.TsSize Es.TsNum])

% Initilize
ind  = 1;
time = 0;
history = [];
testres = [];
% Main loop
while (FlagStop==0) && (time<Es.TimeMax)
	Vnext = Ps.IntegFunc(Vs,Ps,Es);		% Integrate in time
	Vboth = cat(3,Vs,Vnext);			% Concatanate the old and new state, to compare
	[FlagStop,score]=CheckSS(Vboth,Ps,Es);	% Compare the new and old states
    
    if(time<Es.TimeMin) % Don't allow stopping before min-time
        FlagStop=0; 
    end;
    
    if(~isempty(Es.TestFunc)) % make a test on new state if relevant
        testres=PerformTests(Vnext,Ps,Es,Es.TestFunc);
    end;
   
    % update the score (of how much has changed) and the history 
    score = [time log10(score/Es.TimeDst) testres(:)'];
    history = [history; score];
   
    % Move us along
    Vs = Vnext;		
    ind = ind+1;
    time= time + Es.TimeDst;
    
    % Used for Online-drawing, if relevant
    if Es.OlDraw   
           cla;
           subplot(1,2,1);
           Es.PlotFunc(Vs,Ps,Es);
           %plotst(Vs,Ps,Es);
           title(sprintf('step = %d',ind)); 
           subplot(1,2,2);
           plot(history(:,1),history(:,Es.OlDraw+1));
           title(testtext); 
           drawnow; 
    end;
end;
if(time>=Es.TimeMax && ~Es.NoWarning)
    warning(sprintf('Did not converge at time = %d.',time));
end;

% fix up the history we'll return
sz=size(history,2);
if(sz>2)  % Don't save the convergence info if other info exists
    history(:,2)=[];
else      % Or bring it back from the log scale
    history(:,2)=10.^history(:,2);
end;
VsOut = Vs;



%%% SUB FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% stop, pause and finish buttons for OnLine draw option %
   function stopb(~,~)
        FlagStop=1;
   end
   function pauseb(~,~)
        pause();
   end
   function finishb(~,~)
        Es.OlDraw=0;
   end
end


