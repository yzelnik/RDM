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

% For each value to check if/when we reached
for ii=1:length(Es.ReachVal)
    % normalize
    if(initval>Es.ReachVal(1))
        tmpval=-Es.ReachVal(ii);
        tmphist=-bfhist(:,2);
    else
        tmpval=Es.ReachVal(ii);
        tmphist=bfhist(:,2);
    end;
    
    if(initval==Es.ReachVal(1))
        reachtime(ii)=0;    % We started from this value
    elseif(max(tmphist)<tmpval)
        reachtime(ii)=inf;  % Value is never reached
    else
        %[~,minind]=min(abs(bfhist(:,2)-Es.ReachVal(1)));
        crossloc=find(abs(diff(sign(tmphist-tmpval)))>eps); % first sign switch
        reachtime(ii)=bfhist(crossloc(1),1);
    end;

end;

end




