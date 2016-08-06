function [frames,history]=runframes(Vs,Ps,Es,varargin)
% Run integrator (to time Es.Tdest) and get (Es.Frames) snapshots of the system
% [frames,history]=runframes(Vs,Ps,Es)
% Alternatively, if Es.Frames is a vector, then is specifies the time point of the snapshots
% if Es.DynPrm exists, than it allows running different frames with different parameters (Ps). 
% Es.DynPrm is a cell array of these parameters names, while Es.DynVal has the values (in either an array or cell array).
% Es.FramesChoice can be used to specify which frames to get (default is all), with either indices or a logical array format.
% Ps.IntegFunc is the integrator function to use

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Vs = Vs(:,:,1); % Take only first state, if there are more than one.

if((~isfield(Es,'Tdest')) && (length(Es.Frames)>1))
	Es.Tdest = max(Es.Frames);
end;

if(~isfield(Es,'OLdraw'))
    Es.OLdraw=0;
end;

% Deal with dynamic parameters if necessary
if(~isfield(Es,'DynPrm'))
	Es.DynPrm = [];
end;
% Make into cell array if it is not empty
if((~iscell(Es.DynPrm)) && ~isempty(Es.DynPrm))
	Es.DynPrm = {Es.DynPrm};
end;
% Check for a recurring function (operates every-so-often)
if(~isfield(Es,'RecurFunc'))	
    Es.RecurFunc = [];
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


% Calculate any matrices and other auxiliary data before run
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);

% Initilize
if (length(Es.Frames)>1)	% Deal with Es.Frames being a vector
	num = length(Es.Frames);
	Es.Frames = [0 ; Es.Frames(:)];
	Es.Tdest = Es.Frames(end);
	jumps = Es.Frames(2:end)-Es.Frames(1:end-1); 

    if(min(jumps)<0)
       jumps = Es.Frames(2:end);
        warning('Es.Frames should be given with absolute values for time, not relative. Values given were converted.');
    end;
else
	num = Es.Frames;
	jumps = ones(1,num)*Es.Tdest/(num+0.0);
end;

% Es.FramesChoice Allows only some frames to be taken out. Default is to return all. 
% Es.FramesChoice can be given as [1 1 0 0 1 0] or [1 2 5] for getting the first, second and fifth out of six frames
if(~isfield(Es,'FramesChoice') || ((length(Es.FramesChoice)==1) && (Es.FramesChoice(1)==0)))
    Es.FramesChoice = ones(num,1);
elseif (length(Es.FramesChoice) < num)
	temp = zeros(num,1);
	temp(Es.FramesChoice) = 1;
	Es.FramesChoice = temp;
end;

% Setup a recurring function information (function that operates every-so-often)
% Es.RecurFunc can be given as [1 1 0 0 1 0] or [1 2 5] for getting the first, second and fifth out of six frames

if(~isempty(Es.RecurFunc))
    if(~isfield(Es,'RecurFrames') || ((length(Es.RecurFrames)==1) && (Es.RecurFrames(1)==0)))
        Es.RecurFrames = zeros(num,1); % Turn off RecurFunc
    elseif (length(Es.RecurFrames) < num) || (max(Es.RecurFrames)>1)
        temp = zeros(num,1); % Setup specific frames to use RecurFunc
        temp(Es.RecurFrames) = 1;
        Es.RecurFrames = temp;
    end;
end;
%Es.RecurFrames
% Set up GUI for online drawing
FlagStop = 0;
if Es.OLdraw
    clf;
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.6 0.08,0.1],'string','Stop','callback',{@stopb});
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.4 0.08,0.1],'string','Pause','callback',{@pauseb});
    if(isempty(Test))  % If no test was defined, but we are using Online-draw, make something up...
        Test=@T_L2Norm;
    end;
end


frames = zeros([size(Vs) length(nonzeros(Es.FramesChoice))]);
history = [];

frmind = 1;
 %Es.RecurFrames
 
Vnow = Vs;
% Go over each frame
for index=1:num
    % Deal with dynamic parameters if necessary
	if(~isempty(Es.DynPrm))	
		for indprm=1:length(Es. DynPrm)
			if(iscell(Es.DynVal))
				tempval = Es.DynVal{index,indprm};
			else
				tempval = Es.DynVal(index,indprm);
			end;
			if(ischar(Es.DynPrm{indprm}))
				Ps.(Es.DynPrm{indprm})=tempval;
			else
				Ps.Ds(Es.DynPrm{indprm})=tempval;
			end;
		end;
	end;
	Es.Tdest = jumps(index);
    
    % Run a recurring function, if the time is right  
    if(~isempty(Es.RecurFunc))&&(Es.RecurFrames(index))   
        Vnext = Es.RecurFunc(Vnow,Ps,Es);
        Vnow = Vnext;
    end;
    
	% Run integrator to next frame
	Vnext = Ps.IntegFunc(Vnow,Ps,Es);
	Vnow = Vnext;
    
    % Save a frame for later on, if the time is right  
	if(Es.FramesChoice(index))
		frames(:,:,frmind)=Vnow;
		frmind = frmind+1;
	end;
    
    % Save run history (bf data)
    if(~isempty(Test))  
        testres = Test(Vnext,Ps,Es);
        history = [history;  Es.Frames(index) testres(:)'];
    end;
    
    % Option for breaking the loop
    if(FlagStop~=0)
        break;
    end;
    % Used for Online-drawing, if relevant
    if Es.OLdraw   
           cla;
           subplot(1,2,1);
           plotst(Vnow,Ps,Es);
           titletext = sprintf('step = %d',index);
           if(~isempty(Es.DynPrm) && ~iscell(Es.DynVal))
               titletext = sprintf('%s, %s = %.4f',titletext,Es.DynPrm{1},Es.DynVal(index,1));
           end;
           title(titletext); 
           subplot(1,2,2);
           plot(history(:,1),history(:,Es.OLdraw+1));
           title('log10 of score'); 
           drawnow; 
    end;
end;

%%% SUB FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% stop and pause buttons for OnLine draw option %
   function stopb(~,~)
        FlagStop=1;
   end

   function pauseb(~,~)
        pause();
   end

end
   