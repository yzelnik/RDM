function [BfData]=rungrid(Vs,Ps,Es,varargin)
% Run multiple scenarios in parallel, with differnet parameters
% Use a test function (Es.TestFunc) to get a measure/norm along the branch
% points in parameters Es.BFpar, with values specificed by Es.BFrange

% Update online if necessary, but not the state - only parameters
[~,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
Es.InitActive=0; % Allow states to be updated if necessary

% Check if output should be continously written out to file
WriteFlag = 0;
if(isfield(Es,'BFout') & Es.BFout)
	WriteFlag = 1;
end;

if(~isfield(Es,'RunFunc') || isempty(Es.RunFunc))
    Es.RunFunc = @runflow;      % By default use runflow for each scenario
end;

parname = Es.BFpar;
parrange = Es.BFrange;


BfData=[];
for ii=1:length(parrange)
    %disp(ii)
    %disp([Es.BFpar '=' num2str(Es.BFrange(ii))]);
	% Update paramater
	Ps.(parname) = parrange(ii);
    
	% Run the system to SS
	[~,res] = Es.RunFunc(Vs,Ps,Es);
	% Analysze this SS
	%res = Es.TestFunc(VsOut,Ps,Es);
	% Add res to the overall results
	BfData = [BfData; Ps.(parname) res(:)'];
	%disp(BfData)
	if(WriteFlag)
		dlmwrite(Es.BFout,BfData);
	end;
end;
