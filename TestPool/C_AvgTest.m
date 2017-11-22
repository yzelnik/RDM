function [avgt,stdt] = C_AvgTest(Input,Ps,Es,varargin)
% Average out the output of a test function (and calculate it's std as well)

% Update online if necessary
if(nargin>3) [~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'TestFields',[],'CalcRange',[]);

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

if(isempty(Es.TestFields))
 	Es.TestFields=2:size(bfdata,2); % default is all fields but the first one
end;

% What range of frames are we doing the calculation on?
if(isempty(Es.CalcRange)) % by default, range is all data
    calcrange = 1:size(bfdata,1);
else
    if(Es.CalcRange(1)<1)  % is range given in relative size?
        Es.CalcRange = ceil(Es.CalcRange*size(bfdata,1));
    end;
    % is only the first point given? (second point is the end by default)
    if(length(Es.CalcRange)<2) 
        Es.CalcRange(2) = size(bfdata,1);
    end;
    % define range accordingly
    calcrange = Es.CalcRange(1):Es.CalcRange(2);
 
end;

%calcualte the avg & std
avgt=mean(bfdata(calcrange,Es.TestFields),1);
stdt=std(bfdata(calcrange,Es.TestFields),1);

end




