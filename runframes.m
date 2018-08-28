function [frames,history]=runframes(Vs,Ps,Es,varargin)
% Run integrator (to time Es.TimeDst) and get (Es.Frames) snapshots of the system
% [frames,history]=runframes(Vs,Ps,Es)
% Alternatively, if Es.Frames is a vector, then is specifies the time point of the snapshots
% if Es.DynPrm exists, than it allows running different frames with different parameters. 
% Es.DynPrm is a cell array of these parameters names, while Es.DynVal has the values (in either an array or cell array).
% Es.FramesChoice can be used to specify which frames to get (default is all), with either indices or a logical array format.
% Es.RecurFunc and Es.RecurFrames can be used to run a modifying function after before of the frames
% If Es.RecurFrames has more than one column, each column corresponds to a function in Es.RecurFrames


% Default first extra input is for the frames to run
if(~mod(nargin,2)) varargin = ['Es.Frames' varargin]; end;
    
% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'DynPrm',[],'RecurFunc',[],'OlDraw',0,'TestFunc',[],'TsMode','none','FileOut',[],'BfPrm',[],'PlotFunc',@plotst,'RepDynVal',0);
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

% Calculate time step automatically if relevant
[Vs,Ps,Es]=SetupTimeStep(Vs,Ps,Es);
% Calculate any matrices and other auxiliary data before run
[Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);

Vs = Vs(:,:,1); % Take only first state, if there are more than one.

% Initilize
if (length(Es.Frames)>1)	% Deal with Es.Frames being a vector
	num = length(Es.Frames);
	Es.Frames = Es.Frames(:);
	Es.TimeDst = Es.Frames(end);
else
    if(~(Es.Frames==round(Es.Frames)))
        error('If Es.Frames is a single number, it must be an integer.');
    end;
	num = Es.Frames;
    Es.Frames = (1:num)'*Es.TimeDst/num;
    
end;
Es.Frames=round(Es.Frames/Es.TsSize)*Es.TsSize;
jumps = round(diff([0;Es.Frames])/Es.TsSize)*Es.TsSize; 
if(min(jumps)<0)
    error('Es.Frames should be an ordered list')
end;

% Es.FramesChoice allows only some frames to be taken out. Default is to return all. 
% Es.FramesChoice can be given as [1 1 0 0 1 0] or [1 2 5] for getting the first, second and fifth out of six frames
if(~isfield(Es,'FramesChoice') || ((length(Es.FramesChoice)==1) && (Es.FramesChoice(1)==0)))
    Es.FramesChoice = ones(num,1);
elseif (length(Es.FramesChoice) < num) || (max(Es.FramesChoice)>1)
	temp = zeros(num,1);
	temp(Es.FramesChoice) = 1;
	Es.FramesChoice = temp;
end;

% Wrap in cell array if needed
if(~isempty(Es.RecurFunc)&&~iscell(Es.RecurFunc))  
    Es.RecurFunc={Es.RecurFunc};
end;

% Setup a recurring function information (function that operates every-so-often)
% Es.RecurFrames can be given as [1 1 0 0 1 0] or [1 2 5] for getting the first, second and fifth out of six frames
if(~isempty(Es.RecurFunc))
    recfunctype=zeros(length(Es.RecurFunc),1);
    % Check type of recurring functions
    for ii=1:length(Es.RecurFunc)
        txt = func2str(Es.RecurFunc{ii});
        if(strcmp(txt(1:2),'U_'))
            recfunctype(ii)=1;
        end;
    end;
    if(~isfield(Es,'RecurFrames'))
        Es.RecurFunc = []; % Turn off RecurFunc as there are no frames
        Es.RecurFrames = [];
    elseif(size(Es.RecurFrames,1)>1 && size(Es.RecurFrames,2)>1) % if there's more than one set of recur func/frames
        if(~(size(Es.RecurFrames,2)==length(Es.RecurFunc)))
            error('If column number of Es.RecurFrames is more than 1, it should match number of funcs in Es.RecurFunc');
        end;
        Es.RecurFramesExtra=Es.RecurFrames; % Keep the real information here
        Es.RecurFrames = -sum(Es.RecurFrames,2); % just to know if this frames has any recuring func
    else 
        if(length(Es.RecurFrames)==1)
            if(Es.RecurFrames(1)==0)
                % Turn off RecurFunc as there are no frames
                Es.RecurFunc = []; 
            else % assume that this sole value is every how many frames to run RecurFunc 
                temp = zeros(num,1); % Setup specific frames to use RecurFunc
                temp(Es.RecurFrames:Es.RecurFrames:length(Es.Frames)) = 1;
                Es.RecurFrames = temp;
            end;
        elseif(length(Es.RecurFrames)==num-1)
            Es.RecurFrames=[0;Es.RecurFrames]; % assume this is a binary yes/no values per frame
        elseif(length(Es.RecurFrames) < num) || (max(Es.RecurFrames)>1)
            temp = zeros(num,1); % Setup specific frames to use RecurFunc
            temp(Es.RecurFrames) = 1;
            Es.RecurFrames = temp;
        end;
    end;
end;

% If no test was defined, but we are using online-draw, or history was
% requested, make something up...
if(isempty(Es.TestFunc)) && (Es.OlDraw || nargout>1)
    Es.TestFunc=@T_L2Norm;
end;

% Set up GUI for online drawing
FlagStop = 0;
if Es.OlDraw
    clf;
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.7 0.09,0.1],'string','Stop','callback',{@stopb});
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.5 0.09,0.1],'string','Pause','callback',{@pauseb});
    uicontrol('style','pushbutton','units','norm','position',[0.01 0.3 0.09,0.1],'string','Finish','callback',{@finishb});
    testtext=regexprep(func2str(Es.TestFunc),'_','-');
end

% Initilize
frames = zeros([size(Vs) length(nonzeros(Es.FramesChoice))]);
if(~isempty(Es.TestFunc))  
    history = [Es.Frames(:) zeros(size(Es.Frames(:)))];
end;
frmind = 1;
Vnow = Vs;

if(Es.RepDynVal) % replicate Es.DynVal as necessary? (periodic forcing etc...)
    dynmax = size(Es.DynVal,1);
else
    dynmax = num; % no repeat (default)
end;

% Main loop
for index=1:num
    % Deal with dynamic parameters if necessary
	if(~isempty(Es.DynPrm))	
        dynindex = mod(index-1,dynmax)+1;
        % Allow for simple vector change if only 1 parameter is used
        if(length(Es.DynPrm)==1) && (size(Es.DynVal,2)>1)
            tempvals = {Es.DynVal(dynindex,:)}; 
        else
            tempvals = Es.DynVal(dynindex,:);
        end;
        [Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,tempvals,Es.DynPrm,length(Es.BfPrm));
	end;
    
	Es.TimeDst = jumps(index); 
    % Run a recurring function, if the time is right  
    if(~isempty(Es.RecurFunc))&&(Es.RecurFrames(index)) 
        if(Es.RecurFrames(index)>0) % In each frame, run all funcs together or none at all
            for funcind=1:length(Es.RecurFunc)
                if(recfunctype(funcind)==0)
                    Vnext = Es.RecurFunc{funcind}(Vnow,Ps,Es);
                else
                    [Vnext,Ps,Es] = Es.RecurFunc{funcind}(Vnow,Ps,Es);
                end;
                Vnow = Vnext;
            end;
        else % Each function has a differnet timing set
            for funcind=1:length(Es.RecurFunc)
                if(Es.RecurFramesExtra(index,funcind))
                    if(recfunctype(funcind)==0)
                        Vnext = Es.RecurFunc{funcind}(Vnow,Ps,Es);
                    else
                        [Vnext,Ps,Es] = Es.RecurFunc{funcind}(Vnow,Ps,Es);
                    end;
                    Vnow = Vnext;
                end;
            end;
        end;
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
        %history = [history;  Es.Frames(index+1) testres(:)'];
        history(index,2:(length(testres)+1))=testres(:)';
    end;
    
    % Used for Online-drawing, if relevant
    if Es.OlDraw  
           cla;
           if(~isempty(Es.TestFunc))
                subplot(1,2,2);
                plot(history(1:index,1),history(1:index,Es.OlDraw+1));
                title(testtext); 
                subplot(1,2,1);
           end;
           Es.PlotFunc(Vnow,Ps,Es);
           titletext = sprintf('step = %d',index);
           if(~isempty(Es.DynPrm) && ~iscell(Es.DynVal))
               titletext = sprintf('%s, %s = %.4f',titletext,Es.DynPrm{1},Es.DynVal(dynindex,1));
           end;
           title(titletext); 
           drawnow; 
    end;
    if(~isempty(Es.FileOut))
        save(Es.FileOut);
    end;
    % Option for breaking the loop
    if(FlagStop~=0)
        frames = frames(:,:,1:frmind-1); % chop off the end of the data since it is blank
        if(~isempty(Es.TestFunc))
            history = history(1:frmind-1,:);
        end;
        break;
    end;
end;

% fix-up frames
if((frmind-1)<size(frames,3))
    frames=frames(:,:,1:frmind-1);
end;

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
   