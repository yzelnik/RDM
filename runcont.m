function [VsOut,BfData]=runcont(Vs,Ps,Es,varargin)
% Follow points along a branch 
% [BfData,VsOut]=runcont(Vs,Ps,Es)
% Use some steady-state finding function (Es.SsFunc), with NewtonLoop as default
% Use a test function (Es.TestFunc) to get a measure/norm along the branch
% Branch is followed along parameter Es.BfPrm, in points specificed by Es.BfRange

if(~mod(nargin,2)) error('No default extra-input exists for runcont.'); end;
    
% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'BfMaxDiff',0,'TestFunc',@T_L2Norm,'ContUpdate',[],'ContMin',0,'Verbose',0);

% By default, a stopping val is checked on the first column
if(length(Es.BfMaxDiff)<2) ||  (Es.BfMaxDiff(2)<2)
    Es.BfMaxDiff(2)=1;
end;
% In case we want to run several test functions
if(iscell(Es.TestFunc))  
    Es.TestList=Es.TestFunc;
    Es.TestFunc=@T_MultiTest;
end;
% Check if output should be continously written out to file
WriteFlag = 0;
if(isfield(Es,'BfOut') & Es.BfOut)
	WriteFlag = 1;
end;

DefBifPoints = 100;

% get rid of grid-cells for Es.BfPrm
if(iscell(Es.BfPrm)) 
    parname=Es.BfPrm{1};
    % Assume we want to use the other BfPrm's for further analysis inside each run
    Es.BfPrm=Es.BfPrm(2:end);    
else
    parname = Es.BfPrm;
    Es.BfPrm= {};
end;


% get rid of grid-cells for Es.BfRange
if(iscell(Es.BfRange))  
    if(length(Es.BfRange)>1)
        error('Mult-parameter continuation (%s , %s) not supported',Es.BfRange{1},Es.BfRange{2});
    else
        Es.BfRange=Es.BfRange{1};
    end;
end;

if(size(Es.BfRange,2)==1)
    Es.BfRange=Es.BfRange';
end;
if(size(Es.BfRange,2)>4)
    parrange = Es.BfRange(:);
elseif (size(Es.BfRange,2)<3) % assume we only got first and last value
    if(isfield(Es,'BfSmall') & Es.BfSmall)
        parrange = (Es.BfRange(1):Es.BfSmall*((diff(Es.BfRange)>0)*2-1):Es.BfRange(2))';
    else % Or use some arbitary default value
        parrange = (Es.BfRange(1):(Es.BfRange(2)-Es.BfRange(1))/(DefBifPoints-1):Es.BfRange(2))';
    end;
else % Or, we got first&last point, and total-number of points
    parrange = (Es.BfRange(1):(Es.BfRange(2)-Es.BfRange(1))/(Es.BfRange(3)-1):Es.BfRange(2))';
end;

if(~isfield(Es,'SsFunc') || isempty(Es.SsFunc))
	Es.SsFunc = @runnewt;
	Es.OlDraw = 0;
end;


BfData=[];
ii=1;
stopflag=0;

ExtData=[];

while((ii<=length(parrange)) && (~stopflag))
    Vs = Vs + rand(size(Vs))*Es.StSmall*1e-4;
	% Update paramater
    tmpval = parrange(ii);
    if(isempty(strfind(parname,'.')) && isempty(strfind(parname,'('))) 
        Ps.(parname) = tmpval;  % for parameters in Ps (Prefereable)
    else
        eval(sprintf('%s=tmpval;',parname));   % Use eval, less safe/stable
	end;
    if(~isempty(Es.ContUpdate)) % allow for updates in each cont step
        [Vs,Ps,Es]=Es.ContUpdate(Vs,Ps,Es);
    end;
    
%	Ps.(parname) = parrange(ii);
 
    % Run the system to SS
	[VsOut,ExtData] = Es.SsFunc(Vs,Ps,Es);
    if(isempty(VsOut)) % if the result is the bif-data itself
        BfData = [BfData; tmpval ExtData(:)'];
    else
    	% Analysze this SS
        res = Es.TestFunc(VsOut,Ps,Es);
        % Add res to the overall results
        BfData = [BfData; tmpval res(:)'];
    end;
    %disp(BfData(end,:));
    
    if(Es.BfMaxDiff(1) && (ii>1)) % is the new state outside of our preferred bounds?
        if(abs(diff(BfData(ii-1:ii,Es.BfMaxDiff(2)+1)))>Es.BfMaxDiff(1))
            stopflag=1;
        end;
    end;
    % If there's a minimal "success" value, check if we passed the threshold
    if(Es.ContMin) && (ExtData>Es.ContMin)
        stopflag=1;
    end;
	% Update from the last result
	if(stopflag)
        BfData(end,:)=[];
        VsOut=Vs; % push us one state back
    else
        if(~isempty(VsOut)) % update the new state if it's not empty
            Vs = VsOut;
        else
            if(~isempty(Es.BfPrm)) % update the other parameters, if relevant
                [Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,ExtData(1:length(Es.BfPrm)));
            end;
        end;
    end;
	if(WriteFlag)
		dlmwrite(Es.BfOut,BfData,'precision',5);
	end;
    if(Es.Verbose)
        disp(sprintf('step %d: %s=%.4f, ext: %f',ii,parname,parrange(ii),ExtData(1)));
    end;
    ii=ii+1;
	
end;
