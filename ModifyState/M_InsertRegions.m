function VsOut=M_InsertRegions(Vs,Ps,Es,varargin)
% Insert a number of regions in a given state

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if ~isfield(Es,'RegsParms')  
    Es.RegsParms = 1;       % By default, add one (positive) region
end;

Es.RegsParms = [Es.RegsParms(:)' 0 0];  % buffer with zeros
% Vind -> Main variable to work on
if isfield(Es,'Vind')                  
    Es.Vind = Es.Vind(1);
else
    Es.Vind = 1;
end;

% What are the min and max values of all of this?
[md,mn,mx]=T_MinMax(Vs,Ps,Es);      

% How many peaks, and are they positive or negative?
pospeak=1;
regnum=Es.RegsParms(1);
if(regnum<0)
    pospeak=0;
    regnum=-regnum;
end;

% Set value of deleted peaks to?
if(Es.RegsParms(2)==0)
    if(pospeak)
        defval = mn(Es.Vind);
    else
        defval = mx(Es.Vind);
    end;
else
    defval = Es.RegsParms(2);
end;


st = Vs(:,Es.Vind);
if(~pospeak)    % Negate there (and in the end) to simplify reg search
    st=-st;
    defval=-defval;
end;

stdel=M_DeleteRegions(st,Ps,Es,'Es.RegsParms',[1 0]); % Get state minus 1 peak
stdif=M_CenterSt(st-stdel,Ps,Es);         % Single region (that will be added)


defmin = min(Vs(:,Es.Vind))+md(Es.Vind)*1e-3;
nnsm=NeighborSM(1,Ps,Es);

for ii=1:regnum         % Go over regions
    % Init
    %finreg=false(size(st,1),1);
    %bwreg=finreg;
    %plot(min(st,defmax));
    %pause;
    tmpst = max(st,defmin);
    [~,loc]=min(tmpst);    % Find max (starting) point
    %loc
    %st(loc)
    %plotst(st,Ps,Es);
    %pause;
    st = st+M_ShiftSt(stdif,Ps,Es,'Es.ShiftParms',[-loc+Ps.Nx/2 0]); 
    %bwreg(loc)=1;
    
    %while(sum(bwreg)>0)
    %    nn=logical((1-finreg).*(nnsm*bwreg));               % Find nearest neighbors
    %    newreg=logical(nn.*(tmpst<=max(tmpst(logical(bwreg)))));  % those that are smaller go in
    %    finreg(bwreg)=1;                                    % update the final region array
    %    bwreg=newreg;                                       % switch old for the new
        %plot([bwreg(51:100) finreg(51:100)]);
        
        %disp([max(st(logical(bwreg))) nonzeros(nn.*st)'])
        %pause;
        %[newreg' ; bwreg'; finreg']
    %end;
    %st(finreg)=defval;       % delete region
end;

if(~pospeak)    % Negate there (and before loop) to simplify reg search
    st=-st;
end;

VsOut=Vs;
VsOut(:,Es.Vind)=st;

end
