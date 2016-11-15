function regs=SegmentRegions(Vs,Ps,Es,varargin)
% Segment the state into a number of regions
% Threshold is either given by Es.SegThresh, or if not, 
% it is chosen as the average value between min and max values in variable
% Regions above threshold are given ascending integer values (1,2,3...)
% While regions below threhsold are given descending integer values (-1,-2,-3...)
% Values of 0 imply the function did not work properly

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

vmax = max(Vs(:,Es.VarInd(1),1));
vmin = min(Vs(:,Es.VarInd(1),1));
if(isfield(Es,'SegThresh'))
   thresh = Es.SegThresh;
else
    % Find max and min values of the relevant variable in the state, and define the threshold
    thresh = (vmax+vmin)/2;
end;

regs=zeros(size(Vs));
if((vmax-vmin)>Es.StSmall)
    % Use threshold to seperate the state into positive and negative parts (segmentation into bw)
    
    pos=(Vs(:,Es.VarInd(1),1)>thresh);
    neg=(Vs(:,Es.VarInd(1),1)<thresh);
    % define nearest neighbor matrix
    nn=NeighborSM(1,Ps,Es);	

    % Find locations, choose a submatrix of nearest neighbors for positive only, and then find its components
    fpos=find(pos);
    nnpos=nn(fpos,fpos);
    comppos=components(nnpos);

    % Find locations, choose a submatrix of nearest neighbors for negative only, and then find its components
    fneg=find(neg);
    nnneg=nn(fneg,fneg);
    compneg=components(nnneg);

    % Combine the 2 sets of components found
    regs(fpos)=comppos;
    regs(fneg)=-compneg;
end;

end


% Side note: The main part after the definition of pos and neg, can be replaced by:
%   pos=(reshape(Vs(:,Es.VarInd(1),1),Ps.Ny,Ps.Nx)-vmin>thresh);
%   neg=(vmax-reshape(Vs(:,Es.VarInd(1),1),Ps.Ny,Ps.Nx)>thresh);
%   regpos=bwlabel(pos,4);
%   regneg=bwlabel(neg,4);
%   regs=reshape(regpos-regneg,Ps.Ny*Ps.Nx,1);
% Which is generally much faster (X10 or so), but not general
% For example, it does not include periodic boundary conditions
