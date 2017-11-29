function Es=SortOutBfParameters(Es)
% Deal with Es.BfPrm and Es.BfRange, putting it in standard format
% This function is called by runpar
% Three formats are supported for Es.BfRange, for prm # of N (in Es.BfPrm):
% 1) Specific points. N columns, each row is a point in parameter space,
%    each column is for a different parameter. Column min length is 6.
% 2) Array of size 3 or 4. Format is: [LowVal, HighVal, NumVal, Type]
%    Where NumVal gives the number of evaluations, and Type is either:
%    Type = 0: regular grid spacing (default value), 
%    Type < 0: uniform rand to the power of (-Type) (hence, -1 : uniform-rand)
%    0<Type<1: gaussian-rand with std=2/Type (so Type=1 will include ~95%)
%    Type >=2: logarithmic grid spacing (not random). mult-spacing ~Type/NumVal
%    Type=NaN: Same parition with same randomization as last parameter
%    If NumVal=0, then no new points are formed for this parameter
% 3) Cell array, each cell with format as 1) or 2) above

% Do we have parameters who's values are given in cell-arrays?
Es=InsertDefaultValues(Es,'PrmInCell',zeros(1,length(Es.BfPrm)));

% Wrap in cell-array format for convenience 
if(~iscell(Es.BfPrm))  
    Es.BfPrm={Es.BfPrm};
end;

% Put into cell array for convenience, and make sure there's no empty cells
if(~iscell(Es.BfRange))
    Es.BfRange={Es.BfRange};
end;
emptycells = find(cellfun(@isempty,Es.BfRange)); % find empty cells
Es.BfRange(emptycells) = [];                     % remove these cells

curnum=1; % number of different values of parameter (size of total prm-set)

for ii=1:length(Es.BfRange)
	% If we got a cell-array within cell-array, it means cell-data in inside Es.BfRange
	if(iscell(Es.BfRange{ii})) 
        if(~isfield(Es,'CellData'))
            Es.CellData={}; % Create the Es.CellData if needed
        end;
        curcel = length(Es.CellData)+1;
        curlen = length(Es.BfRange{ii});
        Es.CellData{curcel}=Es.BfRange{ii}; % Move cell-arr into its proper place
        Es.BfRange{ii}=[1 curlen -curcel 0]; % Now make the proper arangements in Es.BfRange
	end;
	if(length(Es.BfRange{ii})>5) % type 1 - just take columns of points as-is 
        if(size(Es.BfRange{ii},1)==1) % force into a column (for 1 row of data)
            Es.BfRange{ii}=Es.BfRange{ii}';
        end;
        parvals = Es.BfRange{ii};
        curnum  = size(Es.BfRange{ii},1);
        parnum  = size(Es.BfRange{ii},2);
        tmp     = [0 0 curnum];
    else % type 2
        parnum  = 1; % only one parameter set is specified per cell
        tmp = [Es.BfRange{ii}(:)' 0]; % Pad with zero (for regular grid)
        
        if(tmp(3)<0) % This a "special" type of parameter, given in celldata (-tmp(3))
            Es.PrmInCell(ii)=-tmp(3); % "pointer" to where the data is located
            curnum=tmp(2);
        elseif(tmp(3)>0) 
            curnum=tmp(3); % number of new points multipled by tmp(3)
        end;
        
        % Now determine the type of values for this parameter (grid, rand, etc.)
        if(tmp(4)==0) % regular grid spacing
            parvals = (tmp(1):(tmp(2)-tmp(1))/(curnum-1):tmp(2))';
        elseif(tmp(4)<0) % in general unif^-(tmp(4)), if tmp(4)==-1 then uniform random distribution
            parvals = rand(curnum,1).^(-tmp(4))*(tmp(2)-tmp(1))+tmp(1);
        elseif(tmp(4)<1)   % gaussian/normal random distribution for: 0<tmp(4)<1
            % cuttoff at range edges, tmp(4)=0.5 is four std inside range
            parvals = ((randn(curnum,1)*tmp(4)+2)/4*(tmp(2)-tmp(1))+tmp(1));
            parvals = min(max(parvals,min(tmp(1:2))),max(tmp(1:2)));
        elseif(tmp(4)>=2)  % logaritmic spacing (non-random), consecutive multplication ~ tmp(4)/(tmp(3)-1)
            % create series from 0 to 1, each jump larger than the last
            zerotoone=tmp(4).^(((0:(curnum-1))./(curnum-1)))/(tmp(4)-1)-1/(tmp(4)-1);
            parvals = zerotoone'*(tmp(2)-tmp(1))+tmp(1);
        elseif(isnan(tmp(4))) % each point follows seperation from previous parameter (same distances) 
            tmpvals = (parvals-Es.BfRange{ii-1}(1))/(Es.BfRange{ii-1}(2)-Es.BfRange{ii-1}(1)); % normalize
            parvals = tmpvals*(tmp(2)-tmp(1))+tmp(1);
        else
            error('Unidentified type of value distribution for param: "%s" given by %.2f.',Es.BfPrm{ii},tmp(4));
        end;  
        
	end;  
    
    % Add the new set of parameter values to the final result
    if(ii==1)
        totvals=parvals; % For first parameter we need a new array in any case
    else
        if(abs(tmp(3))>0) % From second parameter on, either increase number of points, or keep the same number
            tmpvals = repmat(totvals,curnum,1); % Make new points, new number is old one multiplied by curnum
            %totvals = [tmpvals reshape(repmat(parvals,1,size(totvals,1))',curnum*size(totvals,1),parnum)];
            totvals = [tmpvals reshape(repmat(reshape(parvals,1,curnum*parnum),size(totvals,1),1),curnum*size(totvals,1),parnum)];
        else
            len = size(totvals,1)/curnum;         % No new points, just spread the new prm-set across current set
            totvals = [totvals reshape(repmat(parvals,len,1),curnum*len,1)];
        end;     
    end;
end;

Es.BfVal = totvals;
end
