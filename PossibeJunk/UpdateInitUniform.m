function VsOut = UpdateInitUniform(Vs,Ps,Es,varargin)
% Initiate a state using initial values given in Vs, and a mask given by Es.StMask
% VsOut = InitStateByMask(Vs,Ps,Es,varargin)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Check if we should run an ODE integration until Steady-State, per initial value set
if(isfield(Es,'InitByODE') && (length(Es.InitByODE)>0))
	if(length(Es.InitByODE)==1)
		Es.InitByODE = Es.InitByODE * ones(size(Vs,1),1);
	end;
else
	Es.InitByODE = zeros(size(Vs,1),1);
end;

% Go over each set of initial values, and run the ODE integration
for ind = 1:length(Es.InitByODE)
	if(Es.InitByODE(ind))
		final = FindODESS(Vs(ind,:),Ps,Es);
		VsOut(ind,:) = final;
	else
		VsOut(ind,:) = Vs(ind,:);
	end;
end;

disp(Vs(:)')	% delete?

end
