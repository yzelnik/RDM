function regsize=T_LowReg(Vs,Ps,Es,varargin)
% Returns the size of the region below threshold Es.SegThresh
% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

regsize = mean(Vs(:,Es.VarInd(1),1)<Es.SegThresh(1));

end
