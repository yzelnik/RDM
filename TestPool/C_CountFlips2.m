function [fliptimes,flipcnt] = C_CountFlips2(Input,Ps,Es,varargin)
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
windowval = Es.FlipThresh(2)/tjmp;  % normalize the smoothing coefficient

% segment the time-series data
%segdata = smooth(data,windowval)>Es.FlipThresh(1);
initdata = diff(data>Es.FlipThresh(1));
tmpflips = [1; find(initdata) ;length(initdata)];
%find(diff(tmpflips)>windowval)'
actflips = tmpflips(find(diff(tmpflips)>windowval));
%[abs(diff(initdata(actflips)))' ; initdata(actflips(1:end-1))']
if(isempty(actflips))
    actflips2=[];
else
    actflips2 = actflips(logical([1 ;abs(diff(initdata(actflips)))]));
end;
%bfhist(actflips,Es.BfFields(1))
%plot(bfhist(:,Es.BfFields(1)), data,'b',bfhist(tmpflips,Es.BfFields(1)),zeros(size(tmpflips))+0.1,'*g',bfhist(actflips,Es.BfFields(1)),zeros(size(actflips))+0.2,'*r',bfhist(actflips2,Es.BfFields(1)),zeros(size(actflips2))+0.3,'*m')
% how many flips between up and down did we have?
flipcnt = length(actflips2(2:end));
%flipcnt = sum(abs(nonzeros(diff(segdata))));
% what is the frequency?
fliptimes = diff(bfhist(actflips2,Es.BfFields(1)));
%flipfreq = 1./[mean(fliptimes(1:2:end)) mean(fliptimes(2:2:end))];
%if(length(fliptimes)<3)
%    fliptimes=rand(randi(10),1)*1000;
%end;
%plot(bfhist(:,Es.BfFields(1)), data,'b', bfhist(:,Es.BfFields(1)),segdata,'r');


end
