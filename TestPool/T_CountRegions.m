function regnum=T_CountRegions(Vs,Ps,Es,varargin)
% Returns the number of different regions in a given state
% Region segmentation is done on the variable given by Es.VarInd
% A threshold of the average of max and min values of the variable is used

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if (~isfield(Es,'VarInd'))
    Es.VarInd = 1;
end;

% Segment the state into positive and negative regions
regs   = SegmentRegions(Vs,Ps,Es);

% The max and min values correspond to the number of positive and negative regions
posnum = max(regs(:));
negnum = -min(regs(:));

% Return these
regnum = [posnum negnum];

end
