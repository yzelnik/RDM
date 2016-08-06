function regsize=T_ZeroStateSize(Vs,Ps,Es,varargin)
% Returns the size of the (largest) zero-domain in a given state 
% The size is given as the fraction of the whole system
% Region segmentation is done on the variable given by Es.Vind
% A threshold used is the geometric mean of Es.STsmall and max(abs) value of state

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;


if (~isfield(Es,'Vind'))
    Es.Vind = 1;
end;


thresh = geomean([Es.STsmall max(abs(Vs(:,Es.Vind)))]);
% Segment the state into positive and negative regions
regs = SegmentRegions(abs(Vs),Ps,Es,'Es.SegThresh',thresh);


% Get Area count for zero-dominated domain
temp = regionprops(-regs.*(regs<0),'Area');

zerosizes=cat(1, temp.Area);

% Return the largest of each, divided by the system size
if(length(temp)) % Make sure there's data here
    regsize = [max(zerosizes)]/(Ps.Nx*Ps.Ny);
else
    regsize = [NaN];
end;

% Alternative (simpler) version, just return precentage that's close to 0
%regsize = mean(abs(Vs(:,Es.Vind))<thresh);

end
