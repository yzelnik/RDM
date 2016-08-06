function [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin)
% Updates parameters online (while running a given function)
% [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin)
% It is assumed the first input is the name of the field in text form
% And the second is its value. For example:
% SomeFunc(...,'Ps.Lx',100)

% Initialization
setnum=floor((nargin-2)/2);
MinSysLen = 10;

% For each couple of input values
for ii=1:setnum

	% Read field name and seperate it
	[Var,Field]=strread(varargin{ii*2-1},'%[^.].%s');
	% Read the value of this field
	Value = varargin{ii*2};
	
	% Give warnings if the name does not make sense. Otherwise update
	if(isempty(Field))
		warning(sprintf('UpdateParameters: Warning. No Field given for %s.',Var{1}))
	else
		if(strcmpi(Var,'Ps'))
			Ps.(Field{1})=Value;
		else 
			if(strcmpi(Var,'Es'))
				% NEED CHANGE - ADD the possibilty of subfield (also for Ps, above)
				Es.(Field{1})=Value;
			else
				warning(sprintf('UpdateParameters: Warning. No Variable by the name of %s.',Var{1}))
			end;
		end;	
	end;
end;
%Es
% Deal with partial or missing Vs
if((~isfield(Es,'InitActive') || Es.InitActive==0) && (size(Vs,1)<MinSysLen) && (size(Vs,1)>0) && (~iscell(Vs)) && ~isempty(Vs))
	Es.InitActive = 1;	% Mark flag so we don't go into an infinite loop
	if(~isfield(Es,'InitFunc'))	
		Es.InitFunc = @M_InitRndSt;	% Default function
	end;
	if(size(Vs,2)<Ps.Vnum)	% Replicate values, for the lazy amongst us
		tmp = repmat(Vs,1,Ps.Vnum);
		Vs = tmp(:,1:Ps.Vnum);
	end;
	Vs = Es.InitFunc(Vs,Ps,Es);
end;


end
