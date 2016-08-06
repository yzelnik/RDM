function [frames,history]=GetSnapshots(Vs,Ps,Es,IntFunc,varargin)
% Run integrator (to time Es.Tdest) and get (Es.Frames) snapshots of the system
% [frames,history]=GetSnapshots(Vs,Ps,Es,IntFunc)
% Alternatively, if Es.Frames is a vector, then is specifies the time point of the snapshots
% if Es.DynPrm exists, than it allows running different frames with different parameters (Ps). 
% Es.DynPrm is a cell array of these parameters names, while Es.DynVal has the values (in either an array or cell array).
% Es.FramesChoice can be used to specify which frames to get (default is all), with either indices or a logical array format.
% IntFunc is the integrator function to use

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if((~isfield(Es,'Tdest')) && (length(Es.Frames)>1))
	Es.Tdest = max(Es.Frames);
end;
% Deal with dynamic parameters if necessary
if(~isfield(Es,'DynPrm'))
	Es.DynPrm = [];
end;
% Make into cell array if it is not empty
if((~iscell(Es.DynPrm)) && length(Es.DynPrm))
	Es.DynPrm = {Es.DynPrm};
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


frames = zeros([size(Vs) length(nonzeros(Es.FramesChoice))]);
history = [];

frmind = 1;
 
Vnow = Vs;
% Go over each frame
for index=1:num
	if(length(Es.DynPrm))	% Deal with dynamic parameters if necessary
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
	% Run integrator to next frame
	Vnext = IntFunc(Vnow,Ps,Es);
	Vnow = Vnext;
	if(Es.FramesChoice(index))
		frames(:,:,frmind)=Vnow;
		frmind = frmind+1;
	end;
    if(~isempty(Test))  % Save run history
        testres=Test(Vnext,Ps,Es);
        history = [history; testres(:)'];
    end;
end;
