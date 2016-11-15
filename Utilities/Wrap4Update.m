function out=Wrap4Update(Vs,Ps,Es,RunFunc,varargin)
% Wrap any function (RunFunc) and run it after UpdateParameters
% out=Wrap4Update(Vs,Ps,Es,RunFunc,varargin)

% Update online if necessary
if(nargin>4) 
    [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); 
end;

%varargout=RunFunc(Vs,Ps,Es);
out=RunFunc(Vs,Ps,Es);

%if(length(varargout)>1)
%    out = varargout;
%else
%    varargout
%    out = varargout{1};
%end;

end
