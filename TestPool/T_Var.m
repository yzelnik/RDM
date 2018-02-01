function varres = T_Var(Vs,Ps,Es,varargin)
% Calculate the variance of Es.VarInd variable(s) (def of Es.VarInd = 1)
% avgvals = T_AvgVal(Vs,varargin)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if (~isfield(Es,'VarInd'))
    Es.VarInd = 1;
end;

for ii=1:length(Es.VarInd)
   varres(ii) = var(Vs(:,Es.VarInd(ii),1));
end

end