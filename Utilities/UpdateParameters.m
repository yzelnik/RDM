function [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin)
% Updates parameters online (while running a given function)
% [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin)
% It is assumed the first input is the name of the field in text form
% And the second is its value. For example:  SomeFunc(...,'Ps.Lx',100)

if(~isfield(Es,'NoWarning'))
    Es.NoWarning=0;
end;

% Initialization
setnum=floor((nargin-2)/2);

% For each couple of input values
for ii=1:setnum

	% Read field name and seperate it
	[Var,Field]=strread(varargin{ii*2-1},'%[^.].%s');
	% Read the value of this field
	Value = varargin{ii*2};
	
	% Give warnings if the name does not make sense. Otherwise update
	if(isempty(Field) && ~Es.NoWarning)
		warning(sprintf('No Field given for %s.',Var{1}))
	else
		if(strcmpi(Var,'Ps'))
			Ps.(Field{1})=Value;
		else 
			if(strcmpi(Var,'Es'))
				% NEED CHANGE - ADD the possibilty of subfield (also for Ps, above)
				Es.(Field{1})=Value;
            elseif (~Es.NoWarning)
				warning(sprintf('No Variable by the name of %s.',Var{1}))
			end;
		end;	
	end;
end;


end
