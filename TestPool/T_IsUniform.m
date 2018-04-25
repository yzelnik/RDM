function [unfflag,maxdif]=T_IsUniform(Vs,Ps,Es,varargin)
% Check if the state is uniform (use Es.UnfThresh as threhsold, or Es.StSmall as default)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Put in some default values of Es
Es=InsertDefaultValues(Es,'VarInd',1,'UnfThresh',[]);

maxdif = max(Vs(:,Es.VarInd(1),1))-min(Vs(:,Es.VarInd(1),1));

if(isempty(Es.UnfThresh))
    unfflag = (maxdif<Es.StSmall);
else
    unfflag = (maxdif<Es.UnfThresh);
end;

end

