function VsOut=M_DeleteRegions(Vs,Ps,Es,varargin)
% Delete a number of regions in a given state

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if ~isfield(Es,'RegPrm')  
    Es.RegPrm = 1;       % By default, delete one (positive) region
end;

Es.RegPrm = [Es.RegPrm(:)' 0 0];  % buffer with zeros
% VarInd -> Main variable to work on
if isfield(Es,'VarInd')                  
    Es.VarInd = Es.VarInd(1);
else
    Es.VarInd = 1;
end;

% What are the min and max values of all of this?
[md,mn,mx]=T_MinMax(Vs,Ps,Es);      

% How many peaks, and are they positive or negative?
pospeak=1;
regnum=Es.RegPrm(1);
if(regnum<0)
    pospeak=0;
    regnum=-regnum;
end;

% Set value of deleted peaks to?
if(Es.RegPrm(2)==0)
    if(pospeak)
        defval = mn(Es.VarInd);
    else
        defval = mx(Es.VarInd);
    end;
else
    defval = Es.RegPrm(2);
end;

st = Vs(:,Es.VarInd);
if(~pospeak)    % Negate there (and in the end) to simplify reg search
    st=-st;
    defval=-defval;
end;

defmax = max(Vs(:,Es.VarInd))-md(Es.VarInd)*1e-1;
nnsm=NeighborSM(1,Ps,Es);

for ii=1:regnum         % Go over regions
    % Init
    finreg=false(size(st,1),1);
    bwreg=finreg;
    %plot(min(st,defmax));
    %pause;
    tmpst = min(st,defmax);
    [~,loc]=max(tmpst);    % Find max (starting) point
    bwreg(loc)=1;
    
    while(sum(bwreg)>0)
        nn=logical((1-finreg).*(nnsm*bwreg));               % Find nearest neighbors
        newreg=logical(nn.*(tmpst<=max(tmpst(logical(bwreg)))));  % those that are smaller go in
        finreg(bwreg)=1;                                    % update the final region array
        bwreg=newreg;                                       % switch old for the new
      
    end;
    st(finreg)=defval;       % delete region
end;

if(~pospeak)    % Negate there (and before loop) to simplify reg search
    st=-st;
end;

VsOut=Vs;
VsOut(:,Es.VarInd)=st;

end
