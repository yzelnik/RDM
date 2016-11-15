function [tscale,timeused] = C_CalcConvergenceTimescale(Input,Ps,Es,varargin)
% Estimate the time scale of convergence (assuming exponential decay)
% tscale=C_CalcConvergenceTimescale(Input,Ps,Es)
% Es.BfFields chooses the [time,norm] fields (def = [1,2])
% input can be either a  bif-list, or set of states (using E.TestFunc)
% tscale returns two time-scale estimates, fit to exp & time-to-half
% tused is the amount of time used for the first tscale calculation


% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'BfFields'))
 	Es.BfFields=[1,2];
end;

if(size(Input,3)>1)     % If we get state data, run the test on each state to form a history
    bfhist(1:size(Input,3)) = 1:size(Input,3);
    for ii=1:size(Input,3)
        temp = Es.TestFunc(Input(:,:,ii),Ps,Es);
        bfhist(ii,2:1+length(temp))=temp;
    end;
else                    % Or, assume we got a history
    bfhist = Input;
end;

time = bfhist(:,Es.BfFields(1));
% Assume the end result is what we converged to
data = abs(bfhist(end,Es.BfFields(2))-bfhist(:,Es.BfFields(2)));

% Find the max difference over time, assume convergence starts after it
[~,ind1] = max(data);
ind1=1;
% find time to half-change from max-change to final value
[~,tmp]  = min(abs(data - (data(ind1)+data(end))/2));
tscale(2)= time(tmp)-time(ind1);

% add-hoc correction, don't use first half of data & ignore the very end
ind1 = max([ind1 round(size(bfhist,1)/2)]);     
ind2 = length(time) - ceil((length(time)-ind1)/10);

% fit to a linear function (with the log of the values in the y axis)
coeffs   = polyfit(time(ind1:ind2),log(data(ind1:ind2)),1);
% use the slope (of a logarithmic function to begin with) to get the timescale and also put in the overall time
tscale(1)= -1/coeffs(1);

timeused = time(ind2)-time(ind1);
% this plot compares the log-of-data with the fit to exp(-t/a)
%plot(time(ind1:ind2),log(data(ind1:ind2)),'r',time(ind1:ind2),coeffs(2)+coeffs(1)*time(ind1:ind2),'k--')
%plot(time(ind1:ind2),(data(ind1:ind2)),'r',time(ind1:ind2),exp(coeffs(2)+coeffs(1)*time(ind1:ind2)),'k--')

end