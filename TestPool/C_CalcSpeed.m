function speed = C_CalcSpeed(Input,Ps,Es,varargin)
% Estimate the speed of change in some test function
% Ad hoc parameters [minpoints,minfaction,smallchange]
% with default [2,10,0.1] can be defined via Es.CalcSpeedPrm

% Update online if necessary
if(nargin>3) [~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'BfFields',[1,2],'CalcSpeedPrm',[2,10,0.1]);

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

%plot(bfhist(:,1),bfhist(:,2))

% ad-hoc constants
minpoints = Es.CalcSpeedPrm(1); % def = 2
minfaction = Es.CalcSpeedPrm(2); % def = 10;
smallchange = Es.CalcSpeedPrm(3); % def = 0.1;

% calculating some min's&max's
buffersize = round(size(bfhist,1)/minfaction);
mintime = max(minpoints,buffersize);
tmp = bfhist(:,Es.BfFields(2));

bad = find(abs(diff(tmp))>smallchange*(max(tmp)-min(tmp)));
for ii=1:length(bad)    % Try to get rid of points where there is a jump thats "too big"
    tmp(bad(ii)+1:end) = tmp(bad(ii)+1:end)-(tmp(bad(ii)+1)-tmp(bad(ii)));
end;
%plot(smooth(tmp(~isnan(tmp)),mintime))

actregion = abs(diff(smooth(tmp(~isnan(tmp)),mintime)));
minact = max(actregion(mintime:end-mintime))*smallchange;
%disp([minact mintime buffersize])

if(isempty(nonzeros(actregion))||isempty(minact))
    begloc = 1; % If no change is noticeable, just give the first and last points
    endloc = length(actregion);
else
    %disp([size(actregion) size(minact)])
    begloc = max([mintime find(actregion>minact, 1, 'first')]);
    endloc = min([size(bfhist,1)-mintime find(actregion>minact, 1, 'last')]);
    endloc = max([begloc+1 endloc]); % Might be redundant?
end;

if((endloc-begloc)>(2*buffersize))
    endloc=endloc-buffersize;   % Try not to use the end if not needed
end;

% calculate the total time and change in values, and with that the speed
if((endloc-begloc)>mintime)
    tottime = bfhist(endloc,Es.BfFields(1))-bfhist(begloc,Es.BfFields(1));
    totval  = bfhist(endloc,Es.BfFields(2))-bfhist(begloc,Es.BfFields(2));
    %[begloc endloc tottime totval]
    %[bfhist(endloc,Es.BfFields(2)) bfhist(begloc,Es.BfFields(2))]
    speed=totval/tottime;
    %speed=(finsz-begsz)/(tlen);
else
    speed=NaN;    % Not enough data
end;

end




