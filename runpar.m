function [StData,BfData]=runpar(Vs,Ps,Es,varargin)
% Run multiple scenarios in parallel, with differnet parameters
% Use a some function (Es.RunFunc) to get a measure/norm (default=runflow)
% Parameters to change are Es.BFpar, with values specificed by Es.BFrange
% Es.BFpar is a string or a cell-array of strings. 
% The string is a parameter name in Ps (e.g. "gamma"), 
% or anything from Ps/Es  with the full hieracry (e.g. "Ps.Ds(2)")
% Three formats are supported for Es.BFrange, for parm # of N (in Es.BFpar):
% 1) Specific points. N columns, each row is a point in parameter space
% 2) N columns of size 3 or 4. Format is: [LowVal; HighVal; NumVal; Type]
%    Where NumVal gives the number of evaluations, and Type is either:
%    Type = 0: regular grid spacing (default value), 
%    Type < 0: uniform rand to the power of (-Type) (hence, -1 : uniform-rand)
%    0<Type<1: gaussian-rand with std=2/Type (so Type=1 will include ~95%)
%    Type > 2: logarithmic grid spacing (not random). mult-spacing=Type/NumVal
%    Type=NaN: Same parition with same randomization as last parameter
%    If NumVal=0, then no new points are formed for this parameter
% 3) Cell array. N cell arrays, each per parameter, with format as 2) above

% Update online if necessary, but not the state - only parameters
if(nargin>3) [~,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es.InitActive  = 0; % Allow states to be updated if necessary
Es.MergeBfData = 1; % If/when using runflow, take all bif data as one row

if(~isfield(Es,'WriteFreq') || isempty(Es.WriteFreq))
    Es.WriteFreq = 100;
end;
if(~isfield(Es,'RunFunc') || isempty(Es.RunFunc))
    Es.RunFunc = @runflow;      % By default use runflow for each scenario
end;
if(~iscell(Es.BFpar))   % Wrap in cell array if not already in one
    Es.BFpar={Es.BFpar};
end;

% Get parameter values sorted out before initiating runs
Es=SortOutBfParameters(Es);

% Choose which (of totvals) parm-value combinations should be run
[partrun,whichruns] = CheckPartialRuns(Es);

% Deal with output to file issues (both .mat and .csv/.txt formats)
[WriteFlag,FileName]= CheckOutput(Es,partrun);

% Initilize before main loop
BfData=[];
StData=[];
writeind = 1;

for ii=whichruns
   	% Update paramaters
    [Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,Es.BFparvals(ii,:));
   
    % Run the system
	[st,bf] = Es.RunFunc(Vs,Ps,Es);
    bf = bf(:)';
    
    % Add data to final result
	BfData(size(BfData,1)+1,1:(length(bf)+size(Es.BFparvals,2))) = [Es.BFparvals(ii,:) bf];
    
    StData = [StData; {st}];
	%disp(BfData)
    
	if(WriteFlag)   % Write to file if needed
        if(writeind>=Es.WriteFreq || ii==whichruns(end))
            if(WriteFlag==1)
                save(FileName,'BfData','StData','Es','Ps');
            else
                dlmwrite(FileName,BfData);
            end;
            writeind=1;
        else
            writeind = writeind+1;  
        end;
	end;
end;

end

%%%%%%%%%%%%%%%%%%%% AUX FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [partrun,whichruns] = CheckPartialRuns(Es)
% Choose which (of Es.BFparvals) parm-value combinations should be run

totvals = Es.BFparvals; % list of parameter-combinations

if(isfield(Es,'RunsChoice') & Es.RunsChoice)
    
	if(length(Es.RunsChoice)>2) % old version: whichruns = Es.RunsChoice;
        whichruns = Es.RunsChoice;  
    else        % Es.RunsChoice = [partition-index partition-total-number]
        whichruns = ceil(size(totvals,1)*(Es.RunsChoice(1)-1)/Es.RunsChoice(2))+1 : ceil(size(totvals,1)*Es.RunsChoice(1)/Es.RunsChoice(2));
    end;
    %disp(whichruns)
    partrun = 1;
else
    whichruns = 1:size(totvals,1); % Default option, for most non-cluster runs
    partrun = 0;
end;

end


function [WriteFlag,FileName]=CheckOutput(Es,partrun)
% Deal with output to file issues (both .mat and .csv/.txt formats)
WriteFlag = 0; % Check if output should be continously written out to file
if(isfield(Es,'FileOut') & Es.FileOut)
	WriteFlag = 1;
    FileEnd='';
    if(isnumeric(Es.FileOut))
        FileName = sprintf('runpar_Of%s_par%s_out%d',func2str(Es.RunFunc),Es.BFpar{1},Es.FileOut);
    else
        if(length(Es.FileOut)>4) % If a predefined file type is used 
            if(strcmp(Es.FileOut(end-3:end),'.mat') || strcmp(Es.FileOut(end-3:end),'.csv') || strcmp(Es.FileOut(end-3:end),'.txt'))
                FileEnd=Es.FileOut(end-3:end);
                FileName=Es.FileOut(1:end-4);
                if(~strcmp(Es.FileOut(end-3:end),'.mat'))
                    WriteFlag=2;    % Only write bif data a text
                end;
            else
                FileName=Es.FileOut;
            end;
        else
            FileName=Es.FileOut;
        end;
    end;
    if(partrun) % If this run is only part of a larger set, change name accordingly
        FileName = sprintf('%s_part%dof%d',FileName,Es.RunsChoice(1),Es.RunsChoice(2));
    end;
    FileName = [FileName FileEnd];
else
    FileName='';
end;

end
