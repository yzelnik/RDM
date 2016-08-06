function [VsOut,history,time]=RunToSS(Vs,Ps,Es,IntFunc,varargin)
% Run integrator until steady state is reached

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

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
	Es.Tstep = EvaluateTS(Vs,Ps,Es,IntFunc)/2;
end;

% Maybe delete? recalculate any matrices and other helping data before run
%Es.fmod=-1; 
%Ps=Ps.SpaFunc([],Ps,Es); 
%Es.fmod=0;

%Es.Tstep = EvaluateTS(Vs,Ps,Es,IntFunc)/2;
%Es.Tstep = 0.5.*(Ps.Lx./Ps.Nx).^2./max(Ps.Ds)./2;
%Tnorm = 1; %0.5.*(Ps.Lx./Ps.Nx).^2./min(Ps.Ds);
Es.Tdest = Es.Tstep*Es.Nsteps;
V = Vs;


%disp(FlagDraw)
FlagStop = 0;
if Es.OLdraw
%    figure;
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.5 0.07,0.1],'string','Stop','callback',{@stopb});
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.3 0.07,0.1],'string','Pause','callback',{@pauseb});
end
%recalculating Es.SSthresh for step size...

Es.SSthresh=Es.SSthresh*Es.Tdest;
%disp([ Es.Tdest Es.Tstep Es.Nsteps])
ind = 1;
time= 0;
history = [];
testres=[];
while (FlagStop==0) && (time<Es.Tmax);
	Vnext = IntFunc(V,Ps,Es);		% Integrate in time
	Vboth = cat(3,V,Vnext);			% Concatanate the old and new state, to compare
	[FlagStop,score]=CheckSS(Vboth,Ps,Es);	% Compare the new and old states
    if(~isempty(Test))
        testres=Test(Vnext,Ps,Es);
    end;
   
    score = [time log10(score/Es.Tdest) testres(:)'];
    %disp(score)
    history = [history; score];
    
    V = Vnext;				
    % Used for Online-drawing
    if Es.OLdraw   
           cla;
           subplot(1,2,1);
           plotst(V,Ps,Es);
           title(sprintf('step = %d',ind)); 
           subplot(1,2,2);
           plot(history(:,1),history(:,Es.OLdraw+1));
           title('log10 of score'); 
           drawnow; 
    end;
    %disp(ind)
    ind = ind+1;
    time= time + Es.Tdest;
end;
if(time>=Es.Tmax)
    warning(sprintf('did not converge at time = %d.',time));
end;

VsOut = V;


%% stop and pause buttons for OnLine draw option
   function stopb(~,~)
        FlagStop=1;
   end

   function pauseb(~,~)
        pause();
   end

end
