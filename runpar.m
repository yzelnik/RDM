function [StData,BfData]=runpar(Vs,Ps,Es,varargin)
% Run multiple scenarios in parallel, with differnet parameters
% Use a some function (Es.RunFunc) to get a measure/norm (default=runflow)
% Parameters to change are Es.BfPrm, with values specificed by Es.BfRange
% Es.BfPrm is a string or a cell-array of strings. 
% The string is a parameter name in Ps (e.g. "gamma"), 
% or anything from Ps/Es  with the full hieracry (e.g. "Ps.Ds(2)")
% Three formats are supported for Es.BfRange, for parm # of N (in Es.BfPrm):
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

if(~mod(nargin,2)) error('No default extra-input exists for runpar.'); end;
    
% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'WriteFreq',100,'RunFunc',@runflow,'RandSeed',0,'FileOut',0,'WriteSt',1);



Es.InitActive  = 0; % Allow states to be updated if necessary
Es.MergeBfData = 1; % If/when using runflow, take all bif data as one row

if(~iscell(Es.BfPrm))   % Wrap in cell array if not already in one
    Es.BfPrm={Es.BfPrm};
end;

% setup randomization issues (relevant if parameters are randomly chosen)
if(Es.RandSeed(1)==0)   
    rng('shuffle');         % Get unique (set by time) seed    
else
    rng(Es.RandSeed(1));	% Randomize with a pre-defined seed
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
    [Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,Es.BfVal(ii,:));
   
    % setup randomization (relevant if randomization is used inside Es.RunFunc)
    if(Es.RandSeed(1)==0)   
        rng('shuffle');         % Get unique (set by time) seed    
    else
        rng(Es.RandSeed(1));	% Randomize with a pre-defined seed
    end;
    % Run the system
	[st,bf] = Es.RunFunc(Vs,Ps,Es);
    bf = bf(:)';
    
    % Add data to final result
	BfData(size(BfData,1)+1,1:(length(bf)+size(Es.BfVal,2))) = [Es.BfVal(ii,:) bf];
    
    StData = [StData; {st}];
	%disp(BfData)
    
	if(WriteFlag)   % Write to file if needed
        if(writeind>=Es.WriteFreq || ii==whichruns(end))
            if(WriteFlag==1)
                save(FileName,'BfData','StData','Es','Ps');
            elseif(WriteFlag==2)
                dlmwrite(FileName,BfData);
            else
                save(FileName,'BfData','Es','Ps');
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
% Choose which (of Es.BfVal) parm-value combinations should be run

totvals = Es.BfVal; % list of parameter-combinations

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
        FileName = sprintf('runpar_Of%s_par%s_out%d',func2str(Es.RunFunc),Es.BfPrm{1},Es.FileOut);
    else
        if(length(Es.FileOut)>4) % If a predefined file type is used 
            if(strcmp(Es.FileOut(end-3:end),'.mat') || strcmp(Es.FileOut(end-3:end),'.csv') || strcmp(Es.FileOut(end-3:end),'.txt'))
                FileEnd=Es.FileOut(end-3:end);
                FileName=Es.FileOut(1:end-4);
                if(~strcmp(Es.FileOut(end-3:end),'.mat'))
                    WriteFlag=2;    % Only write bif data as text
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
    if(~(WriteFlag==2)) && (Es.WriteSt==0)
        WriteFlag=3; % Used to write in mat format, but only the bif data
    end;
else
    FileName='';
end;

end
