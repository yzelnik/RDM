function [times,results]=TimeRuns(Vs,Ps,Es,timestepvec,varargin)
% Time different runs with different time steps

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'RunFunc') || isempty(Es.RunFunc))
Es.RunFunc = @runflow;      % By default use runflow for each scenario
end;

if(~isfield(Es,'BFrange'))
   Es.BFrange=[];
end;
runrep = max(1,length(Es.BFrange));

for ii=1:length(timestepvec)
    for jj=1:runrep
        Es.Tstep=timestepvec(ii);
        if(~isempty(Es.BFrange))
            if(iscell(Es.BFrange))
                Ps.(Es.BFpar) = Es.BFrange{jj};
            else
                Ps.(Es.BFpar) = Es.BFrange(jj);
            end;
        end;
        tic;
        outres=Es.RunFunc(Vs,Ps,Es);
        timetmp=toc;
        times(ii,jj)=timetmp;
        results{ii,jj}=outres;
    end;
end;
   


end

