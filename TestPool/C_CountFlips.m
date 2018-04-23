function [flipcnt,flipfreq] = C_CountFlips(Input,Ps,Es,varargin)
% Count number oftimes the system flips between 2 states
% Use Es.FlipThresh = [val,period] where val is the threshold to determine
% when a flip occurs, and period is the time to smooth over

% Update online if necessary
if(nargin>3) [~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'BfFields',[1,2],'CalcRange',[]);


% If we get state data, run the test on each state to form a history
if(size(Input,3)>1)     
    bfhist(1:size(Input,3),1) = 1:size(Input,3);
    for ii=1:size(Input,3)
        temp = Es.TestFunc(Input(:,:,ii),Ps,Es);
        bfhist(ii,2:1+length(temp))=temp;
    end;
else          % Or, assume we got a history
    bfhist = Input;
end;

data = bfhist(:,Es.BfFields(2));  % get data series
tjmp = mean(diff(bfhist(:,Es.BfFields(1)))); % get avegrage time jump
smoothcoeff = Es.FlipThresh(2)/tjmp;  % normalize the smoothing coefficient

% segment the time-series data
segdata = smooth(data,smoothcoeff)>Es.FlipThresh(1);
% how many flips between up and down did we have?
flipcnt = sum(abs(nonzeros(diff(segdata))));
% what is the frequency?
flipfreq = diff(bfhist([1 end],Es.BfFields(1)))/flipcnt;

%plot(bfhist(:,Es.BfFields(1)), data,'b', bfhist(:,Es.BfFields(1)),segdata,'r');


end
