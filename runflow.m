function [StData,BfData]=runflow(Vs,Ps,Es,varargin)
% Run a flow of some functions consecutively, as specified by Es.FuncList
% [StData,BfData]=runflow(Vs,Ps,Es)
% The function will attempt to automatically identify which functions
% accept which input by name, but this can also be specified by Es.FuncSpec
% Values in Es.FuncSpec(:,1) indicate the folowing input/output:
%   1: function takes a state and returns a state
%   2: function takes a state and returns both state and bif data
%   3: function takes a state and returns bif data
%   4: function takes a set of bif data and returns bif data

% Default first extra input is for the list of functions to run
if(~mod(nargin,2)) varargin = ['Es.FuncList' varargin]; end;
    
% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'MergeBfData',0,'TsMode','none');
% Initilize state if necessary
[Vs,Ps,Es]=InitilizeState(Vs,Ps,Es);


if(~isfield(Es,'FuncList') || isempty(Es.FuncList))
    if(isfield(Es,'TestFunc') && ~isempty(Es.TestFunc))
        if(iscell(Es.TestFunc))  % Allow a lazy-access to MultiTest
            Es.TestList=Es.TestFunc;
            Es.TestFunc=@T_MultiTest;
        end;
        Es.FuncList={@run2ss,Es.TestFunc};  % By default run to steady-state and then run the test function
        Es.TestFunc = [];   % Run test just in the end, not throughout run2ss
    else
       error('No Es.FuncList or Es.TestFunc specificed.');
    end;
end;


% Figure out which function on the list takes which input and output
if(~isfield(Es,'FuncSpec') || isempty(Es.FuncSpec))
	for ii=1:length(Es.FuncList)
        txt = func2str(Es.FuncList{ii});
        if(txt(2)=='_') % "Normal" RDM function
            if (txt(1)=='I' || txt(1)=='M')
                Es.FuncSpec(ii,1)=1;      % Input state, output state
            elseif (txt(1)=='T')
                Es.FuncSpec(ii,1)=3;      % Input state, output bifdata
            elseif (txt(1)=='C')
                Es.FuncSpec(ii,1)=4;      % Input bifdata, output bifdata
            else
                error('Function %s specificed in Es.FuncList is not recognized. Use Es.FuncSpec',txt);
            end;
        elseif strcmp(txt(1:3),'run')
            Es.FuncSpec(ii,1)=2;  % RDM function for running scenario: Input state, output state & bifdata
        else
            error('Function %s specificed in Es.FuncList is not recognized. Use Es.FuncSpec',txt);
        end;
    end;
end;
if(Es.FuncSpec(1)==4)   % Make sure we don't start with a bif-only function
    error('Function %s cannot be used as the first function in the flow, as it accepts bifurcation data.',func2str(Es.FuncList{1}));
end;


if((size(Es.FuncSpec,1)==1) && (length(Es.FuncList)>1))  
    Es.FuncSpec = Es.FuncSpec'; % fit into column in necessary 
end;
if(size(Es.FuncSpec,2)<2)   % Need to add in default values of what information to save
    if(sum(Es.FuncSpec(:,1)==4))
        Es.FuncSpec(:,2)=(Es.FuncSpec(:,1)==4);     % If we have Comparison functions, save only their output and tests coming after
       
        lastcalc = find(Es.FuncSpec(:,1)==4,1,'last');
        if(lastcalc < find(Es.FuncSpec(:,1)==3,1,'last'))
            Es.FuncSpec(find(Es.FuncSpec(lastcalc+1:end,1)==3)+lastcalc,2) = 1;
        end;
    else
        if(sum(Es.FuncSpec(:,1)==3))
            Es.FuncSpec(:,2)=(Es.FuncSpec(:,1)==3); % Save output of Test functions
        else
            Es.FuncSpec(:,2)=(Es.FuncSpec(:,1)==2); % Save output of mixed functions
        end;
        Es.FuncSpec(end,2)=1; % Save output from last function in the flow 
    end;
end;

BfData = [];
StData = [];
for ii=1:length(Es.FuncList)    % Go over functions in the flow
  %  disp([ii size(Vs)])
  %  plotst(Vs,Ps,Es);
  % pause;
    if(Es.FuncSpec(ii,1)==0)        % General functions, mostly for updates
        [Vs,Ps,Es] = Es.FuncList{ii}(Vs,Ps,Es);
    elseif(Es.FuncSpec(ii,1)==1)    % Input state, output state
        StOut = Es.FuncList{ii}(Vs,Ps,Es);
        if(Es.FuncSpec(ii,2)>0)
            StData = cat(3,StData,StOut);
        end;
        Vs=StOut;
    elseif(Es.FuncSpec(ii,1)==2)    %  Input state, output state & bifdata
        [StOut,BfOut] = Es.FuncList{ii}(Vs,Ps,Es);
        if(Es.FuncSpec(ii,2)>0)
            StData = cat(3,StData,StOut);
            BfData = CollectBfData(BfData,BfOut,Es.MergeBfData,ii);
            %BfData(size(BfData,1)+(1:size(BfOut,1)),1:1+size(BfOut,2)) = [repmat(ii,size(BfOut,1),1) BfOut];
        end;
        Vs = StOut(:,:,end);
    elseif(Es.FuncSpec(ii,1)==3)     % Input state, output bifdata
        BfOut=[];
        for stind = 1:size(Vs,3)    % In case there are multiple states
            if(Es.FuncSpec(ii,2)<2)
                newbf = Es.FuncList{ii}(Vs(:,:,stind),Ps,Es);
            elseif (Es.FuncSpec(ii,2)<3)
                [tmp1,tmp2] = Es.FuncList{ii}(Vs(:,:,stind),Ps,Es);
                newbf = [tmp1(:)' tmp2(:)'];
            else
                [tmp1,tmp2,tmp3] = Es.FuncList{ii}(Vs(:,:,stind),Ps,Es);
                newbf = [tmp1(:)' tmp2(:)'  tmp3(:)'];
            end;
            BfOut = [BfOut newbf];
        end;
        if(Es.FuncSpec(ii,2)>0)
            BfData = CollectBfData(BfData,BfOut,Es.MergeBfData,ii);
            %BfData(size(BfData,1)+(1:size(BfOut,1)),1:1+size(BfOut,2)) = [repmat(ii,size(BfOut,1),1) BfOut];
        end;
    elseif(Es.FuncSpec(ii,1)==4)    % Input bifdata, output bifdata
        if(Es.FuncSpec(ii,2)<2)
        	BfOut = Es.FuncList{ii}(BfOut,Ps,Es);
        elseif (Es.FuncSpec(ii,2)<3)
            [tmp1,tmp2] = Es.FuncList{ii}(BfOut,Ps,Es);
            BfOut = [tmp1(:)' tmp2(:)'];
        else
            [tmp1,tmp2,tmp3] = Es.FuncList{ii}(BfOut,Ps,Es);
            BfOut = [tmp1(:)' tmp2(:)'  tmp3(:)'];
        end;
        %BfOut = Es.FuncList{ii}(BfOut,Ps,Es);
        if(Es.FuncSpec(ii,2)>0)
             BfData = CollectBfData(BfData,BfOut,Es.MergeBfData,ii);
%            BfData(size(BfData,1)+(1:size(BfOut,1)),1:1+size(BfOut,2)) = [repmat(ii,size(BfOut,1),1) BfOut];
        end;
    end;    
end;

if((~Es.MergeBfData) && (sum((Es.FuncSpec(:,1)>1).*(Es.FuncSpec(:,2)>0))==1))
    BfData = BfData(:,2:end);	% If only one function was used for creating bif data, don't add extra column
end;
% If we did not read anything into StData (due to the values of Es.FuncSpec),
% then save the last state of the system, if anyone wants to see it
if(isempty(StData))
	StData = Vs;
end;
end

%%%%%%%%%%%% AUX FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BfRes=CollectBfData(BfData,BfOut,mergeflag,ind)
    if(mergeflag)   % Just collect all Bf data in one row
        BfRes = [BfData BfOut(:)'];
    else            % Collect BfData in multiple rows
        BfData(size(BfData,1)+(1:size(BfOut,1)),1:1+size(BfOut,2)) = [repmat(ind,size(BfOut,1),1) BfOut];
        BfRes = BfData;
    end;
    
end
