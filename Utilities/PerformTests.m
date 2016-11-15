function tests=PerformTests(Vs,Ps,Es,TestFunc,varargin)
% Run tests on state-data

% Update online if necessary
if(nargin>4) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~iscell(TestFunc))  % Wrap in cell array if relevant
    funcs={TestFunc}; 
else
    funcs=TestFunc;
end;

szs = size(Vs);

% Go over all states, either in cell format, or matrix format
if(iscell(Vs))
    for ii=1:length(Vs(:))
        for jj=1:length(funcs)
            if(~isempty(Vs{ii}))
                tmp = funcs{jj}(Vs{ii},Ps,Es);
                tests(ii,:)=double(tmp(:)');
            end;
        end;
    end;
    tests = reshape(tests,[szs size(tests,2)]);
else
    for ii=1:size(Vs,3)
        for jj=1:length(funcs)
            tmp = funcs{jj}(Vs(:,:,ii),Ps,Es);
            tests(ii,:)=double(tmp(:)');
        end;
    end;
end; 
    

end
