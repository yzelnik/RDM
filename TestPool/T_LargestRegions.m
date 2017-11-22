function regmax=T_LargestRegions(Vs,Ps,Es,varargin)
% Returns the size of the largest regions in a given state 
% The first is of the positive region, the second of the negative one
% The size of both is given as the fraction of the whole system
% Region segmentation is done on the variable given by Es.VarInd
% A threshold of the average of max and min values of the variable is used

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Segment the state into positive and negative regions
regs = SegmentRegions(Vs,Ps,Es);
plot(regs)
minmax = [min(regs) max(regs)];
if(minmax(1)<0)
    negsizes=histcounts(regs,(minmax(1)-0.5):0);
else
    negsizes=[];
end;
if(minmax(2)>0)
    possizes=histcounts(regs,0.5:(minmax(2)+0.5));
else
    possizes=[];
end;
% Get Area count for positive and negative regions
%temp = regionprops(regs.*(regs>0),'Area');
%possizes=cat(1, temp.Area);
%temp = regionprops(-regs.*(regs<0),'Area');
%negsizes=cat(1, temp.Area);

% Return the largest of each, divided by the system size
if(length(negsizes)*length(possizes)) % Make sure there's data here
    regmax = [max(possizes) max(negsizes)]/(Ps.Nx*Ps.Ny);
else
    regmax = [NaN NaN];%[-1 -1];
end;

end
