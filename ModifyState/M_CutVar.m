function VsOut = M_CutVar(Vs,Ps,Es,varargin)
% Cut down a variable to some level, with (randomized) size
% VsOut = M_CutVar(Vs,Ps,Es)
% Using Es.ModPrm = [cutval,vind,sz,loc,wid], where:
% - cutval tells how much of the variable to cut
% - vind is which variable to change (def=0, is to set by Es.VarInd)
% - sz is the relative size of the system to augment (def=0, all of it)
% - loc is the central location of the cut region (def=0)
% - wid is the width of randomization for size to cut (def=0)
%   cutval<0 gives the value (how much) to be cut, 
% 0<cutval<1 gives the proportaion of the state max to cut,
% 1<cutval<2 gives the (cutval-1) proportion of state mean to cut.
%   cutval>2 gives the volume to cut (in pixels).
% if loc<0 then it is decided at random
% if wid>0 it is uniform, else gaussian (|wid| = std)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es.ModPrm = [Es.ModPrm(:)' 0 0 0 0 0];

if(Es.ModPrm(2)==0) % By default change the Es.VarInd variable
    Es.ModPrm(2)=Es.VarInd;
end;
if(Es.ModPrm(3)==0)
    Es.ModPrm(3)=1; % Augment all the system
elseif(Es.ModPrm(3)>1)
    Es.ModPrm(3)=Es.ModPrm(3)/Ps.Lx;
end;

if(Es.ModPrm(4)>=0) % location of cut region?
    edge = Es.ModPrm(4);
else  % random location
    edge = rand(1);
end;

len = Ps.Nx*Ps.Ny;
% how much to cut? 
if(Es.ModPrm(1)>2)  % cut in absolute pixels
    cutval = - Es.ModPrm(1)/len;  
elseif(Es.ModPrm(1)>1)  % cut proporionally to the mean
    cutval = - (Es.ModPrm(1)-1)*mean(Vs(:,Es.ModPrm(2)));
elseif(Es.ModPrm(1)>0) % cut proporionally to the max
    cutval = - Es.ModPrm(1)*max(Vs(:,Es.ModPrm(2)));
else	% cut in absolute terms
    cutval = Es.ModPrm(1);
end;

if(Es.ModPrm(5)==0) % no size randomization
    cutsz=Es.ModPrm(3); 
elseif(Es.ModPrm(5)>0) % uniform distribution
	wid = min([Es.ModPrm(5) Es.ModPrm(3) 1-Es.ModPrm(3)]);
	cutsz = rand(1)*wid*2-wid+Es.ModPrm(3);
else  % gaussian distribution
    wid = min([Es.ModPrm(5) Es.ModPrm(3) 1-Es.ModPrm(3)]);
    cutsz = randn(1)*wid;
    
    %if(abs(cutsz)>-wid) % apply cutoff
    %    cutsz = wid * sign(cutsz);
    %end;
    cutsz = cutsz + Es.ModPrm(3);
    if(cutsz<0) cutsz=0; end;
    if(cutsz>1) cutsz=1; end;
    disp([wid cutsz])
end;

VsOut=Vs;
if(Es.ModPrm(3)==1) % uniform disturbance
    VsOut(:,Es.ModPrm(2))= VsOut(:,Es.ModPrm(2)) + cutval;
else % sz<1 is relative size, sz>1 is obsolute (in "pixels")
    [reg,regsz]=FindLocalRegion([cutsz edge],Ps,Es);
    
    VsOut(reg,Es.ModPrm(2))= VsOut(reg,Es.ModPrm(2)) + cutval/regsz;
end;

% If we know variables are positive, make sure they remain so	
if((isfield(Es,'NonNeg')) && (Es.NonNeg))
	VsOut = max(0,VsOut);
end;

end
