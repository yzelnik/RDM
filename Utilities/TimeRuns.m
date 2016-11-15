function [times,results]=TimeRuns(Vs,Ps,Es,timestepvec,varargin)
% Time different runs with different time steps

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'RunFunc') || isempty(Es.RunFunc))
Es.RunFunc = @runflow;      % By default use runflow for each scenario
end;

if(~isfield(Es,'BfRange'))
   Es.BfRange=[];
end;
runrep = max(1,length(Es.BfRange));

for ii=1:length(timestepvec)
    for jj=1:runrep
        Es.TsSize=timestepvec(ii);
        if(~isempty(Es.BfRange))
            if(iscell(Es.BfRange))
                Ps.(Es.BfPrm) = Es.BfRange{jj};
            else
                Ps.(Es.BfPrm) = Es.BfRange(jj);
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

