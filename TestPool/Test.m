function tests=Test(Vs,Ps,Es,TestFunc,varargin)
% Run integrator until steady state is reached

% Update online if necessary
if(nargin>4) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~iscell(TestFunc)) funcs={TestFunc}; end;
%tests=zeros(1);
szs = size(Vs);

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
