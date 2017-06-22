function reachtime = C_ReachVal(Input,Ps,Es,varargin)
% Calculate the time it takes for a test to reach Es.ReachVal

% Update online if necessary
if(nargin>3) [~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:}); end;

if(~isfield(Es,'BfFields') || isempty(Es.BfFields))
 	Es.BfFields=[1,2];
end;

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

% Make sure that initial value is always lower than Es.ReachVal
initval = bfhist(1,2);
if(initval>Es.ReachVal(1))
    Es.ReachVal=-Es.ReachVal;
    bfhist(:,2)=-bfhist(:,2);
end;

if(initval==Es.ReachVal(1))
    reachtime=0;    % We started from this value
elseif(max(bfhist(:,2))<Es.ReachVal(1))
    reachtime=inf;  % Value is never reached
else
    %[~,minind]=min(abs(bfhist(:,2)-Es.ReachVal(1)));
    crossloc=find(diff(sign(bfhist(:,2)-Es.ReachVal(1))));
    reachtime=bfhist(crossloc(1),1);
end;

end




