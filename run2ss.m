function [VsOut,history,time]=run2ss(Vs,Ps,Es,varargin)
% Run integrator until steady state is reached 
% [VsOut,history,time]=run2ss(Vs,Ps,Es)
% steady-state is set by ES.SSthresh (or, if not defined, by Es.STsmall) 

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
Vs = Vs(:,:,1); % Take only first state, if there are more than one.

IntegFunc = Ps.IntegFunc;

if(~isfield(Es,'Tmax'))
	Es.Tmax = inf;
end;
if(~isfield(Es,'Nsteps'))
	Es.Nsteps = 1e2;
end;
if(~isfield(Es,'SSthresh'))	% Unless otherwise specified, check if in one "time unit" a threshold of STsmall is passed
	Es.SSthresh = Es.STsmall;%*(Es.Tstep*Es.Nsteps);
end;
if(~isfield(Es,'Tmode'))
	Es.Tmode = 'none';
end;
if(~isfield(Es,'OLdraw'))	% Should the results be automatically drawn and updated while integrating?
	Es.OLdraw = 0;
end;
% Deal with test functions (such as L2Norm) that save the history of the run
if(~isfield(Es,'TestFunc'))	
    Test=[];
else
	if(iscell(Es.TestFunc))
        Test=@T_GetStats;
    else
        Test=Es.TestFunc;
    end;
end;

% Calculate time step automatically if relevant
if(strcmp(Es.Tmode,'auto'))
	Es.Tstep = EvaluateTS(Vs,Ps,Es);
end;

% Calculate any matrices and other auxiliary data before run
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);

%Es.Tstep = EvaluateTS(Vs,Ps,Es,IntegFunc)/2;
%Es.Tstep = 0.5.*(Ps.Lx./Ps.Nx).^2./max(Ps.Ds)./2;
%Tnorm = 1; %0.5.*(Ps.Lx./Ps.Nx).^2./min(Ps.Ds);
Es.Tdest = Es.Tstep*Es.Nsteps;


% Set up GUI for online drawing
FlagStop = 0;
if Es.OLdraw
    clf;
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.6 0.08,0.1],'string','Stop','callback',{@stopb});
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.4 0.08,0.1],'string','Pause','callback',{@pauseb});
end

%recalculating Es.SSthresh for step size...
Es.SSthresh=Es.SSthresh*Es.Tdest;
%disp([ Es.Tdest Es.Tstep Es.Nsteps])
ind = 1;
time= 0;
history = [];
testres=[];
while (FlagStop==0) && (time<Es.Tmax);
	Vnext = IntegFunc(Vs,Ps,Es);		% Integrate in time
	Vboth = cat(3,Vs,Vnext);			% Concatanate the old and new state, to compare
	[FlagStop,score]=CheckSS(Vboth,Ps,Es);	% Compare the new and old states
    
    if(~isempty(Test))
        testres=Test(Vnext,Ps,Es);
    end;
   
    score = [time log10(score/Es.Tdest) testres(:)'];
    history = [history; score];
    
    Vs = Vnext;		
    ind = ind+1;
    time= time + Es.Tdest;
    
    % Used for Online-drawing, if relevant
    if Es.OLdraw   
           cla;
           subplot(1,2,1);
           plotst(Vs,Ps,Es);
           title(sprintf('step = %d',ind)); 
           subplot(1,2,2);
           plot(history(:,1),history(:,Es.OLdraw+1));
           title('log10 of score'); 
           drawnow; 
    end;
end;
if(time>=Es.Tmax)
    warning(sprintf('did not converge at time = %d.',time));
end;

sz=size(history,2);
if(sz>2)  % Don't save the convergence info if other info exists
    history(:,2)=[];
else      % Or bring it back from the log scale
    history(:,2)=10.^history(:,2);
end;
VsOut = Vs;



%%% SUB FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% stop and pause buttons for OnLine draw option %
   function stopb(~,~)
        FlagStop=1;
   end

   function pauseb(~,~)
        pause();
   end

end


