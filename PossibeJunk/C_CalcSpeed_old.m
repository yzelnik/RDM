function speed = C_CalcSpeed(Input,Ps,Es,varargin)
% Estimate the speed of change in some test function

% Update online if necessary
[~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});

if(~isfield(Es,'BfFields'))
Es.BfFields=[1,2];
end;

if(size(Input,3)>1)     % If we get state data, run the test on each state to form a history
bfhist(1:size(Input,3),1) = 1:size(Input,3);
for ii=1:size(Input,3)
temp = Es.TestFunc(Input(:,:,ii),Ps,Es);
bfhist(ii,2:1+length(temp))=temp;
end;
else                    % Or, assume we got a history
bfhist = Input;
end;

%size(bfhist)

minpoints = 2;
minfaction = 10;
buffersize = round(size(bfhist,1)/minfaction);
%maxpoints = 1000;

mintime = max(minpoints,buffersize);
actregion = abs(diff(smooth(bfhist(:,Es.BfFields(2)),mintime)));
minact = max(actregion)*1e-2;
%disp([minact mintime buffersize])

%plot(actregion)
%[min(find(actregion)) max(find(actregion))]
begloc = max([mintime find(actregion>minact, 1, 'first')]);
endloc = min([size(bfhist,1)-mintime find(actregion>minact, 1, 'last')]);

if((endloc-begloc)>(2*buffersize))
endloc=endloc-buffersize;
end;
%[begloc endloc]
%if((endloc-begloc)>maxpoints)
%    tmp=(endloc+begloc)/2;
%    begloc = round(tmp-maxpoints/2);
%    endloc = round(tmp+maxpoints/2);
%end;
%plot(actregion>minact)
%find(actregion, 1, 'last')
%[begloc endloc]
%mintime = 5*min(diff(bfhist(:,Es.BfFields(1))));

%finsz = bfhist(end,Es.BfFields(2));
%begsz = bfhist(1,Es.BfFields(2));
%[~,indmax]=max(bfhist(:,Es.BfFields(2)));
%[~,indmin]=min(bfhist(:,Es.BfFields(2)));
%inds = sort([indmin indmax]);
%tlen=max([indmax indmin]);
%tlen = bfhist(inds(2),Es.BfFields(1))-bfhist(inds(1),Es.BfFields(1));
%tval = bfhist(inds(2),Es.BfFields(2))-bfhist(inds(1),Es.BfFields(2));
tottime = bfhist(endloc,Es.BfFields(1))-bfhist(begloc,Es.BfFields(1));
totval  = bfhist(endloc,Es.BfFields(2))-bfhist(begloc,Es.BfFields(2));
%disp([tottime  bfhist(end,Es.BfFields(1))])
%disp([totval tottime endloc begloc])
%plot(bfhist(begloc:endloc,Es.BfFields(1))-bfhist(begloc,Es.BfFields(1)),bfhist(begloc:endloc,Es.BfFields(2)));
%disp([indmin indmax tlen tval])
%[bfhist(endloc,Es.BfFields(2)) bfhist(begloc,Es.BfFields(2))]
%[tottime totval]

%plot(bfhist(:,2:4));
%pause;
if(tottime>mintime)
speed=totval/tottime;
%speed=(finsz-begsz)/(tlen);
else
speed=0;
end;
spd2=speed;
%disp(speed)

minpnts = 4;
maxshare = 0.4; % Max share a point can have in total change, when calculating center-of-mass
dd=diff(bfhist(:,Es.BfFields(2)));
dd2 = dd;
problemdd=find(abs(dd)>(sum(abs(dd))*maxshare));    % points with too high of-a-change
dd2(problemdd)=0;
%plot((dd2(:)'.*(1:length(dd2))))
       % "center-of-mass" of change calculation
       if((~sum(dd2)==0) && (length(nonzeros(dd2))>minpnts))
       half = ceil(minpnts/2);
       %center=round(sum(sqrt(abs(dd(:)')).*(1:length(dd)))/sum(abs(sqrt(dd))))
                                  center=round(sum(dd2(:)'.*(1:length(dd2)))/sum(dd2));
                                                   if((center<half) || (center>length(dd2)+half))
                                                   %warning('center-of-mass calculation failed. using middle.');
                                                   center=round(length(dd2)/2);
                                                   end;
                                                   %sum(dd2)
                                                   
                                                   begloc = center-half;
                                                   endloc = center+half;
                                                   %disp([begloc endloc])
                                                   
                                                   tottime = bfhist(endloc,Es.BfFields(1))-bfhist(begloc,Es.BfFields(1));
                                                   totval  = bfhist(endloc,Es.BfFields(2))-bfhist(begloc,Es.BfFields(2));
                                                   slope(1) = totval/tottime;
                                                   ind=1;
                                                   while((begloc>half) && (endloc<=(length(dd)-half)))
                                                   tottime = bfhist(endloc+half,Es.BfFields(1))-bfhist(begloc,Es.BfFields(1));
                                                   totval  = bfhist(endloc+half,Es.BfFields(2))-bfhist(begloc,Es.BfFields(2));
                                                   slope1 = totval/tottime;
                                                   tottime = bfhist(endloc,Es.BfFields(1))-bfhist(begloc-half,Es.BfFields(1));
                                                   totval  = bfhist(endloc,Es.BfFields(2))-bfhist(begloc-half,Es.BfFields(2));
                                                   slope2 = totval/tottime;
                                                   if(abs(slope(ind)-slope1)<abs(slope(ind)-slope2))
                                                   slope(ind+1)=slope1;
                                                   endloc=endloc+half;
                                                   %disp('right')
                                                   else
                                                   slope(ind+1)=slope2;
                                                   begloc=begloc-half;
                                                   %disp('left')
                                                   end;
                                                   ind=ind+1;
                                                   end;
                                                   speed=slope(end);
                                                   else
                                                   speed=0;
                                                   end;
                                                   
                                                   %disp(speed)
                                                   
                                                   speed = [speed spd2];
                                                   %[begloc endloc]
                                                   %subplot(1,2,1);
                                                   %plot(smooth(slope))
                                                   %subplot(1,2,2);
                                                   %plotbf(bfhist)
                                                   
                                                   %maxdist = min([floor((center-1)/half) floor((length(dd)-center)/half)]);
                                                   %for ii=1:maxdist
                                                   %    begloc = center-ii*half;
                                                   %    endloc = center+ii*half;
                                                   %    tottime = bfhist(endloc,Es.BfFields(1))-bfhist(begloc,Es.BfFields(1));
                                                   %    totval  = bfhist(endloc,Es.BfFields(2))-bfhist(begloc,Es.BfFields(2));
                                                   %    slope(ii) = totval/tottime;
                                                   %end;
                                                   
                                                   end