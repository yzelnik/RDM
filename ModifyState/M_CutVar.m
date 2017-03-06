function VsOut = M_CutVar(Vs,Ps,Es,varargin)
% Cut down a variable to some level
% VsOut = M_CutVar(Vs,Ps,Es)
% Using Es.ModPrm = [cutval,vind,sz,loc], where:
% - cutval tells how much of the variable to cut
% - vind is which variable to change (def=0, is to set by Es.VarInd)
% - sz is the relative size of the system to augment (def=0, all of it)
% - loc is the location from which to start measuring the size (def=0)
% 0<cutval<1 gives the proportaion of the state max to cut,
% or curval<0 just gives how much to cut
% if loc<0 then it is decided at random

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es.ModPrm = [Es.ModPrm(:)' 0 0 0 0];

if(Es.ModPrm(2)==0)
    Es.ModPrm(2)=Es.VarInd;
end;
if(Es.ModPrm(3)==0)
    Es.ModPrm(3)=1; % Augment all the system
elseif(Es.ModPrm(3)>1)
    Es.ModPrm(3)=Es.ModPrm(3)/Ps.Lx;
end;

if(Es.ModPrm(4)>=0)
    edge = Es.ModPrm(4);
else
    edge = rand(1);
end;

len = Ps.Nx*Ps.Ny;
% how much to cut? 
if(Es.ModPrm(1)>0)  % 
    cutval = - Es.ModPrm(1)*max(Vs(:,Es.ModPrm(2)));
else
    cutval = Es.ModPrm(1);
end;

VsOut=Vs;
if(Es.ModPrm(3)>=1)
    VsOut(:,Es.ModPrm(2))= VsOut(:,Es.ModPrm(2)) + cutval;
else
    nnsm=NeighborSM(1,Ps,Es);
    sitenum = ceil(Es.ModPrm(3)*len);
    mean(sum(nnsm));
    % initilize vector of differnet-order neighbors
    tmpvec  = zeros(size(Vs,1),1);
    tmpvec(ceil(edge*len))=1;
    
    counter = 1;
    while counter<sitenum  % iteratively find neighboring sites
        newvec = logical(nnsm*tmpvec);
        newvec = newvec - tmpvec.*newvec;
        tmpvec = tmpvec + newvec;
        counter = sum(tmpvec ~= 0);
        %plot(newvec); pause;
    end;
    
    if(counter>sitenum)  % make sure we have exactly sitenum sites
        inds = find(newvec);
        tmpvec(inds(1:(counter-sitenum)))=0;     
    end;

    locs=logical(tmpvec);
    %locs = mod(ceil(edge*len)+(1:ceil(Es.ModPrm(3)*len))-1,len)+1;
    
    VsOut(locs,Es.ModPrm(2))= VsOut(locs,Es.ModPrm(2)) + cutval/Es.ModPrm(3);
end;

% If we know variables are positive, make sure they remain so	
if((isfield(Es,'NonNeg')) && (Es.NonNeg))
	VsOut = max(0,VsOut);
end;

end
