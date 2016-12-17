function stats = C_MultiCalc(Vs,Ps,Es,varargin)
% Run multiple calc-functions and combine output
% Es.CalcList should have a list of functions to run,

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Put in some default values of Es
Es=InsertDefaultValues(Es,'CalcList',[],'CalcOuts',[]);

if(isempty(Es.CalcList))
    error('Need a Es.CalcList to be defined');
end;
if(~iscell(Es.CalcList)) % Force calc-list into a cell structure
    Es.CalcList={Es.CalcList};
end;

% Allows the calc-outs to be put inside the CalcList.
if(isempty(Es.CalcOuts))  
    % If inside the cell-array there are numbers following the function
    % names, than assume these are the calc-outs
    ind=0; temp={};
    for ii=1:length(Es.CalcList)

        if(isnumeric(Es.CalcList{ii}))
            Es.CalcOuts(ind,1:length(Es.CalcList{ii})) = Es.CalcList{ii}(:)';
        else
            ind = ind+1;  % Assumed to be a function handle
            temp{ind}=Es.CalcList{ii};
        end;
    end;
    Es.CalcList = temp;
end;


num = length(Es.CalcList);
% Make sure there are enough outputs, by putting the default value of 1 per function
if(size(Es.CalcOuts,1)<length(Es.CalcList))
    Es.CalcOuts(length(Es.CalcList),1)=0;
end;
if(size(Es.CalcOuts,2)<2)
    Es.CalcOuts(1,2)=0;
end;
Es.CalcOuts(:,1) = max(1,Es.CalcOuts(:,1));

stats = [];
% Go over different functions
for ii=1:num
	if(Es.CalcOuts(ii,1)==1)
		o1 = Es.CalcList{ii}(Vs,Ps,Es);
        tmp = o1(:)';
	elseif (Es.CalcOuts(ii,1)==2)
		[o1,o2] =  Es.CalcList{ii}(Vs,Ps,Es);
        tmp = [o1(:)' o2(:)'];

	elseif (Es.CalcOuts(ii,1)==3)
		[o1,o2,o3] =  Es.CalcList{ii}(Vs,Ps,Es);
        tmp = [o1(:)' o2(:)' o3(:)'];
		stats = [stats o1(:)' o2(:)' o3(:)'];
	elseif (Es.CalcOuts(ii,1)==4)
		[o1,o2,o3,o4] =  Es.CalcList{ii}(Vs,Ps,Es);
        tmp = [o1(:)' o2(:)' o3(:)' o4(:)'];
		stats = [stats o1(:)' o2(:)' o3(:)' o4(:)'];
	elseif (Es.CalcOuts(ii,1)==5)
		[o1,o2,o3,o4,o5] =  Es.CalcList{ii}(Vs,Ps,Es);
        tmp = [o1(:)' o2(:)' o3(:)' o4(:)' o5(:)'];
		stats = [stats o1(:)' o2(:)' ];
    else
        tmp = [];
		warning('MultiCalc does not support more than 5 outputs per function');
	end;
    if(Es.CalcOuts(ii,2)>0)  % Get only specific outsputs if requested
        tmp = tmp(nonzeros(Es.CalcOuts(ii,2:end)));
    end;
    stats = [stats tmp];
end;

end