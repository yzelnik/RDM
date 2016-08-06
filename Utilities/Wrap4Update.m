function varargout=Wrap4Update(Vs,Ps,Es,RunFunc,varargin)
% Wrap any function (RunFunc) and run it after UpdateParameters
% out=Wrap4Update(Vs,Ps,Es,RunFunc,varargin)

% Update online if necessary
if(nargin>4) 
    [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); 
end;

varargout=RunFunc(Vs,Ps,Es);

end
