function [avgt,stdt] = C_AvgTest(Input,Ps,Es,varargin)
% Average out the output of a test function (and calculate it's std as well)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% If we get state data, run the test on each state to form a bif table
if(size(Input,3)>1)     
    bfdata(1:size(Input,3),1) = 1:size(Input,3);
    for ii=1:size(Input,3)
        temp = Es.TestFunc(Input(:,:,ii),Ps,Es);
        bfdata(ii,2:1+length(temp))=temp;
    end;
else          % Or, assume we got a bif data
    bfdata = Input;
end;

if(~isfield(Es,'TestFields') || isempty(Es.TestFields))
 	Es.TestFields=2:size(bfdata,2); % default is all fields but the first one
end;

%calcualte the avg & std
avgt=mean(bfdata(:,Es.TestFields),1);
stdt=std(bfdata(:,Es.TestFields),1);

end




