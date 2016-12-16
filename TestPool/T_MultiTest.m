function stats = T_MultiTest(Vs,Ps,Es,varargin)
% Run multiple test functions and combine output
% Es.TestList should have a list of functions to run,
% Or otherwise Es.TestFunc could hold such a list
% Es.TestOuts can give the number of outputs to take from each function

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Put in some default values of Es
Es=InsertDefaultValues(Es,'TestFunc',[],'TestList',[],'TestOuts',[]);

if(isempty(Es.TestList))
    if(isempty(Es.TestFunc))
        error('Need a Es.TestList or Es.TestFunc to be defined');
    else
        Es.TestList=Es.TestFunc;
    end;
end;
if(~iscell(Es.TestList)) % Force test-list into a cell structure
    Es.TestList={Es.TestList};
end;

% Allows the test-outs to be put inside the testlist.
if(isempty(Es.TestOuts))  
    % If inside the cell-array there are numbers following the function
    % names, than assume these are the test-outs
    ind=0; temp={};
    for ii=1:length(Es.TestList)

        if(isnumeric(Es.TestList{ii}))
            Es.TestOuts(ind,1:length(Es.TestList{ii}(:))) = Es.TestList{ii}(:)';
        else
            ind = ind+1;  % Assumed to be a function handle
            temp{ind}=Es.TestList{ii};
        end;
    end;
    Es.TestList = temp;
end;

num = length(Es.TestList);
% Make sure there are enough outputs, by putting the default value of 1 per function
if(size(Es.TestOuts,1)<length(Es.TestList))
    Es.TestOuts(length(Es.TestList),1)=0;
end;
if(size(Es.TestOuts,2)<2)
    Es.TestOuts(1,2)=0;
end;
Es.TestOuts(:,1) = max(1,Es.TestOuts(:,1));

stats = [];
% Go over different functions
for ii=1:num
	if(Es.TestOuts(ii,1)==1)
		o1 = Es.TestList{ii}(Vs,Ps,Es);
        tmp = o1(:)';
	elseif (Es.TestOuts(ii,1)==2)
		[o1,o2] =  Es.TestList{ii}(Vs,Ps,Es);
        tmp = [o1(:)' o2(:)'];

	elseif (Es.TestOuts(ii,1)==3)
		[o1,o2,o3] =  Es.TestList{ii}(Vs,Ps,Es);
        tmp = [o1(:)' o2(:)' o3(:)'];
		stats = [stats o1(:)' o2(:)' o3(:)'];
	elseif (Es.TestOuts(ii,1)==4)
		[o1,o2,o3,o4] =  Es.TestList{ii}(Vs,Ps,Es);
        tmp = [o1(:)' o2(:)' o3(:)' o4(:)'];
		stats = [stats o1(:)' o2(:)' o3(:)' o4(:)'];
	elseif (Es.TestOuts(ii,1)==5)
		[o1,o2,o3,o4,o5] =  Es.TestList{ii}(Vs,Ps,Es);
        tmp = [o1(:)' o2(:)' o3(:)' o4(:)' o5(:)'];
		stats = [stats o1(:)' o2(:)' ];
    else
        tmp = [];
		warning('MultiTest does not support more than 5 outputs per function');
	end;
    if(Es.TestOuts(ii,2)>0)  % Get only specific outsputs if requested
        tmp = tmp(nonzeros(Es.TestOuts(ii,2:end)));
    end;
    stats = [stats tmp];
end;

