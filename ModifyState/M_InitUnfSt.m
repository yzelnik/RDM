function VsOut = M_InitUnfSt(Vs,Ps,Es,varargin)
% form a new state using initial values given in Vs
% Es.OdeInit can specify if these values should be run as ODEs to reach steady-state
% VsOut = M_InitUnfSt(Vs,Ps,Es,varargin)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;
% Make sure there's no spatial heterogenous parameter
Ps=ForcePrmMeanField(Ps);

% Should we actually do anything? (is the state ready already?)
if(size(Vs,1)<(Ps.Nx*Ps.Ny))	
	% Check if we should run an ODE integration until Steady-State, per initial value set
	if(isfield(Es,'OdeInit') && ~isempty(Es.OdeInit))
		if(length(Es.OdeInit)==1)
			Es.OdeInit = repmat(Es.OdeInit,size(Vs,1),1);
		end;
	else	% Don't run any integration, just use initial values as given
		Es.OdeInit = zeros(size(Vs,1),1);
	end;
	
	% Go over each set of initial values, and run the ODE integration if necessary
	for ind = 1:length(Es.OdeInit)
		if(Es.OdeInit(ind))
			VsTmp(1,:,ind) = FindODESS(Vs(ind,:),Ps,Es);
		else
			VsTmp(1,:,ind) = Vs(ind,:);
		end;	
	end;
	
	% repeat to form a uniform state
	VsOut = repmat(VsTmp,Ps.Nx*Ps.Ny,1);
else
	VsOut = Vs;
end;

end
