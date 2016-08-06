function [VsOut,BfData]=runcont(Vs,Ps,Es,varargin)
% Follow points along a branch 
% [BfData,VsOut]=runcont(Vs,Ps,Es)
% Use some steady-state finding function (Es.SSfunc), with NewtonLoop as default
% Use a test function (Es.TestFunc) to get a measure/norm along the branch
% Branch is followed along parameter Es.BFpar, in points specificed by Es.BFrange

% Update online if necessary
[~,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

DefBifPoints = 100;

if(isfield(Es,'TestFunc') && iscell(Es.TestFunc))  % Allow a lazy-access to GetStats
	Es.TestList=Es.TestFunc;
	Es.TestFunc=@T_GetStats;
end;


% Check if output should be continously written out to file
WriteFlag = 0;
if(isfield(Es,'BFout') & Es.BFout)
	WriteFlag = 1;
end;
if(~isfield(Es,'BfMaxDiff') || isempty(Es.BfMaxDiff))
	Es.BfMaxDiff = 0;   % max diff to stop continuation
end;

% get rid of grid-cells for Es.BFpar
if(iscell(Es.BFpar)) 
    parname=Es.BFpar{1};
    % Assume we want to use the other BFpar's for further analysis inside each run
    Es.BFpar=Es.BFpar(2:end);    
else
    parname = Es.BFpar;
    Es.BFpar= {};
end;


% get rid of grid-cells for Es.BFrange
if(iscell(Es.BFrange))  
    if(length(Es.BFrange)>1)
        error('Mult-parameter continuation (%s , %s) not supported',Es.BFrange{1},Es.BFrange{2});
    else
        Es.BFrange=Es.BFrange{1};
    end;
end;

if(size(Es.BFrange,2)>4)
    parrange = Es.BFrange(:);
elseif (size(Es.BFrange,2)<3) % assume we only got first and last value
    if(isfield(Es,'BFsmall') & Es.BFsmall)
        parrange = (Es.BFrange(1):Es.BFsmall*((diff(Es.BFrange)>0)*2-1):Es.BFrange(2))';
    else % Or use some arbitary default value
        parrange = (Es.BFrange(1):(Es.BFrange(2)-Es.BFrange(1))/(DefBifPoints-1):Es.BFrange(2))';
    end;
else % Or, we got first&last point, and total-number of points
    parrange = (Es.BFrange(1):(Es.BFrange(2)-Es.BFrange(1))/(Es.BFrange(3)-1):Es.BFrange(2))';
end;

if(~isfield(Es,'SSfunc') || isempty(Es.SSfunc))
	Es.SSfunc = @NewtonLoop;
	Es.OLdraw = 0;
end;


BfData=[];
ii=1;
stopflag=0;

while((ii<=length(parrange)) && (~stopflag))
    %disp(ii);
    %disp([Es.BFpar '=' num2str(Es.BFrange(ii))]);
	Vs = Vs + rand(size(Vs))*Es.STsmall;
	% Update paramater
	Ps.(parname) = parrange(ii);
    
	% Run the system to SS
	[VsOut,ExtData] = Es.SSfunc(Vs,Ps,Es);
    
    if(isempty(VsOut))
        BfData = [BfData; Ps.(parname) ExtData(:)'];
    else
    	% Analysze this SS
        res = Es.TestFunc(VsOut,Ps,Es);
        % Add res to the overall results
        BfData = [BfData; Ps.(parname) res(:)'];
    end;
    %disp(BfData(end,:));
    if(Es.BfMaxDiff && (ii>1))
        if(abs(diff(BfData(ii-1:ii,2)))>Es.BfMaxDiff)
            stopflag=1;
        end;
    end;
    
	% Update from the last result
	if(stopflag)
        BfData(end,:)=[];
    else
        if(~isempty(VsOut))
            Vs = VsOut;
        else
            if(~isempty(Es.BFpar))
                [Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,ExtData(1:length(Es.BFpar)));
            end;
        end;
    end;
	if(WriteFlag)
		dlmwrite(Es.BFout,BfData,'precision',5);
	end;
    ii=ii+1;
end;
