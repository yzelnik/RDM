function frontsz=T_FrontSize(Vs,Ps,Es,varargin)
% Estimate the size of a front

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'VarInd',1,'FrontThresh',0.05);

thresh = Es.FrontThresh; % relative threshold for top and bottom of front
[~,minvals,maxvals] = T_MinMax(Vs);  % get extreme values

for ii=1:length(Es.VarInd) % for each variable asked for
    vind = Es.VarInd(ii);
    % find the location where the threshold are passed
    [~,botind]=min(abs(Vs(:,vind)-(minvals(vind)*thresh+maxvals(vind)*(1-thresh))));
    [~,topind]=min(abs(Vs(:,vind)-(maxvals(vind)*thresh+minvals(vind)*(1-thresh))));
    %[botind topind]
    frontsz(ii)=abs(topind-botind)*Ps.Lx/Ps.Nx; % renormalize to real size
end;


end



