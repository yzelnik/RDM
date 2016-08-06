function stats = T_GetStats(Vs,Ps,Es,varargin)
% Get stat results from several functions combined

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'TestOuts'))
	Es.TestOuts = [];
end;
if(~isfield(Es,'TestList') || isempty(Es.TestList))
    if(~isfield(Es,'TestFunc') || isempty(Es.TestFunc))
        error('Need a Es.TestList or Es.TestFunc to be defined');
    else
        Es.TestList=Es.TestFunc;
    end;
end;
if(~iscell(Es.TestList))
    Es.TestList={Es.TestList};
end;

if(isempty(Es.TestOuts))  % Allow the test-outs to be put inside the testlist.
    ind=0; temp={};
    for ii=1:length(Es.TestList)

        if(isnumeric(Es.TestList{ii}))
            Es.TestOuts(ind,1:length(Es.TestList{ii})) = Es.TestList{ii}(:)';
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
		warning('GetStats does not support more than 5 outputs per function');
	end;
    if(Es.TestOuts(ii,2)>0)  % Get only specific outsputs if requested
        tmp = tmp(nonzeros(Es.TestOuts(ii,2:end)));
    end;
    stats = [stats tmp];
end;

