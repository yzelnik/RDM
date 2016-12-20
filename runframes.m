function [frames,history]=runframes(Vs,Ps,Es,varargin)
% Run integrator (to time Es.TimeDst) and get (Es.Frames) snapshots of the system
% [frames,history]=runframes(Vs,Ps,Es)
% Alternatively, if Es.Frames is a vector, then is specifies the time point of the snapshots
% if Es.DynPrm exists, than it allows running different frames with different parameters (Ps). 
% Es.DynPrm is a cell array of these parameters names, while Es.DynVal has the values (in either an array or cell array).
% Es.FramesChoice can be used to specify which frames to get (default is all), with either indices or a logical array format.

% Default first extra input is for the frames to run
if(~mod(nargin,2)) varargin = ['Es.Frames' varargin]; end;
    
% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'DynPrm',[],'RecurFunc',[],'OlDraw',0,'TestFunc',[],'TsMode','none','FileOut',[]);
% Initilize state if necessary
[Vs,Ps,Es]=InitilizeState(Vs,Ps,Es);

% Allow for specific frames to be defined
if((~isfield(Es,'TimeDst')) && (length(Es.Frames)>1))
	Es.TimeDst = max(Es.Frames);
end;
% Put dynamic parameters into cell array if it is not empty
if((~iscell(Es.DynPrm)) && ~isempty(Es.DynPrm))
	Es.DynPrm = {Es.DynPrm};
end;
% Allow simple-dynamical variable setup (As row instead of column)
if(~isempty(Es.DynPrm) && ~isempty(Es.DynVal) && (size(Es.DynVal,1)==1))
    Es.DynVal=Es.DynVal';
end;
% In case we want to run several test functions
if(iscell(Es.TestFunc))  
    Es.TestList=Es.TestFunc;
    Es.TestFunc=@T_MultiTest;
end;

% Calculate any matrices and other auxiliary data before run
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);

Vs = Vs(:,:,1); % Take only first state, if there are more than one.

% Initilize
if (length(Es.Frames)>1)	% Deal with Es.Frames being a vector
	num = length(Es.Frames);
	Es.Frames = [0 ; Es.Frames(:)];
	Es.TimeDst = Es.Frames(end);
else
	num = Es.Frames;
    Es.Frames = (0:num)*Es.TimeDst/num;
end;
Es.Frames=round(Es.Frames/Es.TsSize)*Es.TsSize;
jumps = round(diff(Es.Frames)/Es.TsSize)*Es.TsSize; 
if(min(jumps)<0)
    error('Es.Frames should be an ordered list')
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
% Es.RecurFrames can be given as [1 1 0 0 1 0] or [1 2 5] for getting the first, second and fifth out of six frames

if(~isempty(Es.RecurFunc))
    if(~isfield(Es,'RecurFrames') || ((length(Es.RecurFrames)==1) && (Es.RecurFrames(1)==0)))
        Es.RecurFrames = zeros(num,1); % Turn off RecurFunc
    elseif (length(Es.RecurFrames) < num) || (max(Es.RecurFrames)>1)
        temp = zeros(num,1); % Setup specific frames to use RecurFunc
        temp(Es.RecurFrames) = 1;
        Es.RecurFrames = temp;
    end;
end;

% Set up GUI for online drawing
FlagStop = 0;
if Es.OlDraw
    clf;
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.6 0.08,0.1],'string','Stop','callback',{@stopb});
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.4 0.08,0.1],'string','Pause','callback',{@pauseb});
    if(isempty(Es.TestFunc))  % If no test was defined, but we are using Online-draw, make something up...
        Es.TestFunc=@T_L2Norm;
    end;
    testtext=regexprep(func2str(Es.TestFunc),'_','-');
end


frames = zeros([size(Vs) length(nonzeros(Es.FramesChoice))]);
history = [];

frmind = 1;

Vnow = Vs;
% Go over each frame
for index=1:num
    % Deal with dynamic parameters if necessary
	if(~isempty(Es.DynPrm))	
        % Allow for simple vector change if only 1 parameter is used
        if(length(Es.DynPrm)==1) && (size(Es.DynVal,2)>1)
            tempvals = {Es.DynVal(index,:)}; 
        else
            tempvals = Es.DynVal(index,:);
        end;
        [Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,tempvals,Es.DynPrm);
        
		%for indprm=1:length(Es.DynPrm)
	    %	 if(iscell(Es.DynVal))  
		%		tempval = Es.DynVal{index,indprm}; % Read from a cell-array
        %    elseif(length(Es.DynPrm)==1)
        %        tempval = Es.DynVal(index,:);      % Read the whole row
        %    else
		%		tempval = Es.DynVal(index,indprm); % Read from a single column
		%	end;
		%	if(ischar(Es.DynPrm{indprm}))
		%		Ps.(Es.DynPrm{indprm})=tempval;
		%	else
		%		Ps.Ds(Es.DynPrm{indprm})=tempval;
		%	end;
		%end;
	end;
	Es.TimeDst = jumps(index);
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
    if(~isempty(Es.TestFunc))  
        testres = PerformTests(Vnext,Ps,Es,Es.TestFunc);
        history = [history;  Es.Frames(index) testres(:)'];
    end;
    
    % Option for breaking the loop
    if(FlagStop~=0)
        break;
    end;
    % Used for Online-drawing, if relevant
    if Es.OlDraw   
           cla;
           subplot(1,2,1);
           plotst(Vnow,Ps,Es);
           titletext = sprintf('step = %d',index);
           if(~isempty(Es.DynPrm) && ~iscell(Es.DynVal))
               titletext = sprintf('%s, %s = %.4f',titletext,Es.DynPrm{1},Es.DynVal(index,1));
           end;
           title(titletext); 
           subplot(1,2,2);
           plot(history(:,1),history(:,Es.OlDraw+1));
           title(testtext); 
           drawnow; 
    end;
    if(~isempty(Es.FileOut))
        save(Es.FileOut);
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
   