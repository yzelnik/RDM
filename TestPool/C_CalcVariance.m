function varout = C_CalcVariance(Input,Ps,Es,varargin)
% Estimate the variance of the result of some test function

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

% What range of frames are we doing the calculation on?
if(isempty(Es.CalcRange)) % by default, range is from middle to end
    calcrange = round(size(bfhist,1)/2):size(bfhist,1);
else
    if(Es.CalcRange(1)<1)  % is range given in relative size?
        Es.CalcRange = ceil(Es.CalcRange*size(bfhist,1));
    end;
    % is only the first point given? (second point is the end by default)
    if(length(Es.CalcRange)<2) 
        Es.CalcRange(2) = size(bfhist,1);
    end;
    % define range accordingly
    calcrange = Es.CalcRange(1):Es.CalcRange(2);
 
end;

% calculate the variance in this range
varout = var(bfhist(calcrange,Es.BfFields(2)));


end
