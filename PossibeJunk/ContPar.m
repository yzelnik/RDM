function [points,VsOut]=ContPar(Vs,Ps,Es,IntFunc,TestFunc,varargin)
% Use an integrator and test function for continuation of a paramterer

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Check if output should be continously written out
WriteFlag = 0;
if(isfield(Es,'BFout') & Es.BFout)
	WriteFlag = 1;
end;
parname = Es.BFpar;
parrange = Es.BFrange;

points=[];
for ii=1:length(parrange)
    %disp(ii);
    disp([Es.BFpar '=' num2str(Es.BFrange(ii))]);
	Vs = Vs + rand(size(Vs))*Es.STsmall;
	% Update paramater
	Ps.(parname) = parrange(ii);
	% Run the system to SS
	VsOut = RunToSS(Vs,Ps,Es,IntFunc);
	% Analysze this SS
	res = TestFunc(VsOut,Ps,Es);
	% Add res to the overall results
	points = [points; Ps.(parname) res(:)'];
	% Use the last result to get faster convergence
	Vs = VsOut;
	
	if(WriteFlag)
		dlmwrite(Es.BFout,points);
	end;
end;
