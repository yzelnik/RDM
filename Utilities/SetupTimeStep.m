function [Vs,Ps,Es]=SetupTimeStep(Vs,Ps,Es,varargin)
% Setup time-stepping
% VsNew=ChangeRes(Vs,Ps,Es,PsNew)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'TsSize',[],'Verbose',0);

if(length(Es.TsMode)<4)
    Es.TsMode=[];
else
    % Calculate time step automatically ?
    if(strcmp(Es.TsMode(1:4),'auto')) 
    	tmp = EvaluateTS(Vs,Ps,Es);
        if(length(Es.TsMode)==4)
            Es.TsSize=tmp;
        elseif(Es.TsMode(5)=='*') % if we want larger ts than auto
            Es.TsSize=tmp*str2num(Es.TsMode(6:end));
        elseif(Es.TsMode(5)=='/') % if we want smaller ts than auto
            Es.TsSize=tmp/str2num(Es.TsMode(6:end));
        else
            error('Automatic time-step should have Es.TsMode be one of {auto/2,auto,auto*1.3}');
        end;
        if(Es.Verbose)
            disp(sprintf('ts = %e',Es.TsSize));
        end;
    elseif(strcmp(Es.TsMode(1:4),'none')) 
        % did someone forget the define the time-step?
        if(isempty(Es.TsSize) || Es.TsSize==0)
            Es.TsSize = EvaluateTS(Vs,Ps,Es);
            if(Es.Verbose)
                disp(sprintf('No timestep was defined. using: ts = %e',Es.TsSize));
            end;
        end;
    else
        error('Value of Es.TsMode is not recognized');
    end;
end;

end

