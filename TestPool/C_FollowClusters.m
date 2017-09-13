function stats = C_FollowClusters(Input,Ps,Es,varargin)
% Follow the birth, merges, growth and death of clusters

% Update online if necessary
if(nargin>3) [~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'BfFields',[1,2]);
small = 1e-8;

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
change = [0; diff(data)];
events = Es.RecurFrames;
if(max(events)>2)
    tmp=zeros(size(Es.Frames));
    tmp(events)=1;
    events = tmp(1:length(Es.Frames))';
end;


birth  = ((change>small)&(events));
growth = ((abs(change)<small)&(events));
merge  = ((change<-small)&(events));
death  = ((change<-small)&(~events));

birth(1)=1;
growth(1)=0;

if(size(bfhist,2)>Es.BfFields(2)) % are we calculating the overall avg?
    avgchange = diff([0;bfhist(:,Es.BfFields(2)+1)]);
    growth = ((avgchange<0)&(growth));
end;

if(size(bfhist,2)>Es.BfFields(2)+1) % are we following the biggest region?
    bigchange = diff([0;bfhist(:,Es.BfFields(2)+2)]);
    bigind  = [0; diff(bigchange)>small];
    bigrowth= (bigind&growth);
    bigmerge= (bigind&merge);
    bigvec  = [sum(bigrowth) sum(bigmerge) mean(bigchange(bigrowth)) mean(bigchange(bigmerge))];
else
    bigvec = [NaN NaN NaN NaN];
end;

%[birth(1) growth(1) merge(1) death(1)]
%plotwf(Input,Ps,Es);

%disp(bigchange(bigmerge))
%hold on;
%plot(repmat([0;1000],1,sum(birth)),repmat(find(birth)'/10,2,1),'r')
%plot(repmat([0;1000],1,sum(growth)),repmat(find(growth)'/10,2,1),'g')
%plot(repmat([0;1000],1,sum(merge)),repmat(find(merge)'/10,2,1),'b')
%plot(repmat([0;1000],1,sum(death)),repmat(find(death)'/10,2,1),'k')
%hold off;

stats = [sum(birth) sum(growth) sum(merge) sum(death) mean(data) max(data) bigvec];

end