function [Vs,Ps,Es]=InitilizeState(Vs,Ps,Es,varargin)
% A Utility function to fill in the Vs input using an ODE solution
% [Vs,Ps,Es]=InitilizeState(Vs,Ps,Es)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

MinSysLen = 10;

% Put in some default values of Es
Es=InsertDefaultValues(Es,'InitActive',0,'InitFunc',@M_InitRndSt);

% Deal with partial or missing Vs
if((Es.InitActive==0) && ~isempty(Vs) && (size(Vs,1)<Ps.Nx*Ps.Ny) && (size(Vs,1)>0) && (~iscell(Vs)) )
	Es.InitActive = 1;	% Mark flag so we don't go into an infinite loop

	if(size(Vs,2)<Ps.VarNum)	% Replicate values, for the lazy amongst us
		tmp = repmat(Vs,1,Ps.VarNum);
		Vs = tmp(:,1:Ps.VarNum);
	end;
	Vs = Es.InitFunc(Vs,Ps,Es);
end;

end


