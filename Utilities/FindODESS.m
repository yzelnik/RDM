function [state,T,Y] = FindODESS(InitVals,Ps,Es,varargin)
% Find an ODE (uniform) Steady-State
% state = FindODESS(InitVals,Ps,Es)

% Update online if necessary
if(nargin>3) [~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:}); end;
% Default values
Es=InsertDefaultValues(Es,'TimeDst',0,'OdeTime',100,'MaxOdeTime',1000,'OdeThresh',1e-6,'NoWarning',0);

% Baseline simulation time
basetime = Es.OdeTime;
 % Make sure we don't have endless loops 
Es.InitActive = 1;     

% Update online if necessary
[~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});

if(isempty(InitVals))
	InitVals = ones(1,Ps.VarNum);
end;
Es.JacMode=0;
%if(Es.NoWarning)
%    opts = odeset('Stats','off');
%else
%    opts = odeset('Stats','on');
%end;

% Go over each set of initial values, and run the ODE integration
for ind = 1:size(InitVals,1)
    start = InitVals(ind,:);
    maxchange = inf;
    % keep going while convergence is not good and we are not at max time
    while(maxchange/basetime>Es.OdeThresh && basetime*2<Es.MaxOdeTime)
        % run ODE integration
        [T,Y] = ode45(@ODE_shell,(0:4)*basetime/4,start,[],Ps,Es);
        % calculate change in last leg of simulation
        if(size(Y,1)<5)
            basetime=inf;
        else
            maxchange = max(abs(Y(4,:)-Y(5,:)));
            basetime=basetime*2;
            start = Y(5,:); % ready for next iteration
        end;
    end;
    if(maxchange/basetime>Es.OdeThresh) && (~Es.NoWarning)
        warning('Solution does not appear to converge.')
	end;
    if(size(Y,1)<5)
        state(ind,:)=NaN;
    else
        state(ind,:) = Y(5,:);
    end;
end;

end



