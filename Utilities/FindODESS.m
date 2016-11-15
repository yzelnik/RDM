function [state,T,Y] = FindODESS(InitVals,Ps,Es,varargin)
% Find an ODE (uniform) Steady-State
% state = FindODESS(InitVals,Ps,Es)

Es=InsertDefaultValues(Es,'TimeDst',0);

slowchange = 0.01;  % Used to test that a solution is converging. 
deftime = max(100,Es.TimeDst);
if((isfield(Es,'OdeTime')) && Es.OdeTime)
    deftime=Es.OdeTime;
end;

Es.InitActive = 1;      % make sure we don't have endless loops 

% Update online if necessary
[~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});

if(isempty(InitVals))
	InitVals = zeros(1,Ps.VarNum);
end;
Es.JacMode=0;

% Go over each set of initial values, and run the ODE integration
for ind = 1:size(InitVals,1)
    % run ODE integration
    [T,Y] = ode45(@ODE_shell,(0:4)*deftime/4,InitVals(ind,:),[],Ps,Es);
    
    if(sum(abs(Y(2,:)-Y(3,:))*slowchange<abs(Y(4,:)-Y(5,:))))
        % It seems we are not converging still.
        [~,Y2] = ode45(@ODE_shell,(0:4)*deftime/4,Y(5,:),[],Ps,Es);
        %disp([abs(Y(1,:)-Y(5,:)) abs(Y2(1,:)-Y2(5,:))])
        if(sum(abs(Y(1,:)-Y(5,:))*slowchange<abs(Y2(1,:)-Y2(5,:))))
            warning('Solution does not appear to converge.')
        end;
    end;
	state(ind,:) = Y(5,:);
end;
%disp(state)
end

