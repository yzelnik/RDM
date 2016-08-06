function VsOut = M_InitUnfSt(Vs,Ps,Es,varargin)
% form a new state using initial values given in Vs
% Es.InitByODE can specify if these values should be run as ODEs to reach steady-state
% VsOut = M_InitUnfSt(Vs,Ps,Es,varargin)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Should we actually do anything? (is the state ready already?)
if(size(Vs,1)<(Ps.Nx*Ps.Ny))	
	% Check if we should run an ODE integration until Steady-State, per initial value set
	if(isfield(Es,'InitByODE') && (length(Es.InitByODE)>0))
		if(length(Es.InitByODE)==1)
			Es.InitByODE = repmat(Es.InitByODE,size(Vs,1),1);
		end;
	else	% Don't run any integration, just use initial values as given
		Es.InitByODE = zeros(size(Vs,1),1);
	end;
	
	% Go over each set of initial values, and run the ODE integration if necessary
	for ind = 1:length(Es.InitByODE)
		if(Es.InitByODE(ind))
			VsTmp(1,:,ind) = FindODESS(Vs(ind,:),Ps,Es);
		else
			VsTmp(1,:,ind) = Vs(ind,:);
		end;	
	end;
	%disp(VsTmp(:)')	% delete?
	
	% repeat to form a uniform state
	VsOut = repmat(VsTmp,Ps.Nx*Ps.Ny,1);
else
	VsOut = Vs;
end;

end
