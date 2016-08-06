function Es=SortOutBfParameters(Es)
% Deal with Es.BFpar and Es.BFrange, putting it in standard format
% This function is called by runpar
% Three formats are supported for Es.BFrange, for parm # of N (in Es.BFpar):
% 1) Specific points. N columns, each row is a point in parameter space
% 2) N columns of size 3 or 4. Format is: [LowVal; HighVal; NumVal; Type]
%    Where NumVal gives the number of evaluations, and Type is either:
%    Type = 0: regular grid spacing (default value), 
%    Type < 0: uniform rand to the power of (-Type) (hence, -1 : uniform-rand)
%    0<Type<1: gaussian-rand with std=2/Type (so Type=1 will include ~95%)
%    Type > 2: logarithmic grid spacing (not random). mult-spacing=Type/NumVal
%    Type=NaN: Same parition with same randomization as last parameter
%    If NumVal=0, then no new points are formed for this parameter
% 3) Cell array. N cell arrays, each per parameter, with format as 2) above


% setup randomization issues (relevant if parameters are randomly chosen)
if(~isfield(Es,'RandSeed'))
    Es.RandSeed=0;
end;
if(Es.RandSeed(1)==0)   
    rng('shuffle');         % Get unique (set by time) seed    
else
    rng(Es.RandSeed(1));	% Randomize with a pre-defined seed
end;

% Wrap in cell-array format for convenience 
if(~iscell(Es.BFpar))  
    Es.BFpar={Es.BFpar};
end;
if(~isfield(Es,'ParInCell')) % Do we have parameters who's values are given in cell-arrays?
    Es.ParInCell=zeros(length(Es.BFpar),1);
end;


% If Es.BFrange is one row (not cell) that's probably just a mistake to fix
if(~iscell(Es.BFrange) && size(Es.BFrange,1)==1)
    Es.BFrange=Es.BFrange';
end;

% If Es.BFrange is in (2) format, change it into (3) format for convenience 
if(~iscell(Es.BFrange) && size(Es.BFrange,1)<5)
    for ii=1:size(Es.BFrange,2)
        tmp{ii}=Es.BFrange(:,ii);
    end;
    Es.BFrange=tmp;
end;

curnum=1; % number of different values of parameter  

if(iscell(Es.BFrange))  % Now if we're in (3) format (or (2), in effect)
    for ii=1:length(Es.BFrange)
        if(iscell(Es.BFrange{ii}))  % if we got a cell-array within cell-array
            if(~isfield(Es,'CellData'))
                Es.CellData={}; % Create the Es.CellData if needed
            end;
            curcel = length(Es.CellData)+1;
            curlen = length(Es.BFrange{ii});
            Es.CellData{curcel}=Es.BFrange{ii}; % move cell-arr into its proper place
            Es.BFrange{ii}=[1;curlen;-curcel;0]; % now make the proper arangements in Es.BFrange.
        end;
        tmp = [Es.BFrange{ii}(:) ; 0]; % Padd with zero (for regular grid)
        
        if(tmp(3)<0) % This a "special" type of parameter, given in celldata (-tmp(3))
            Es.ParInCell(ii)=-tmp(3); % "pointer" to where the data is located
            curnum=tmp(2);
        end
        if(tmp(3)>0) 
            curnum=tmp(3); % number of new points multipled by tmp(3)
        end;
        
        if(tmp(4)==0) % regular grid spacing
            parvals = (tmp(1):(tmp(2)-tmp(1))/(curnum-1):tmp(2))';
        elseif(tmp(4)<0) % in general unif^-(tmp(4)), if tmp(4)==-1 then uniform random distribution
            parvals = rand(curnum,1).^(-tmp(4))*(tmp(2)-tmp(1))+tmp(1);
        elseif(tmp(4)<1)   % gaussian/normal random distribution for: 0<tmp(4)<1
            % cuttoff at range edges, tmp(4)=0.5 is four std inside range
            parvals = ((randn(curnum,1)*tmp(4)+2)/4*(tmp(2)-tmp(1))+tmp(1));
            parvals = min(max(parvals,min(tmp(1:2))),max(tmp(1:2)));
        elseif(tmp(4)>2)    
            % logaritmic spacing (non-random), consecutive multplication = tmp(4)/tmp(3)
            parvals = tmp(4).^(((1:tmp(3))./tmp(3)))/tmp(4)*(tmp(2)-tmp(1))+tmp(1);
        elseif(isnan(tmp(4))) % each point follows seperation from previous parameter (same distances) 
            tmpvals = (parvals-Es.BFrange{ii-1}(1))/(Es.BFrange{ii-1}(2)-Es.BFrange{ii-1}(1)); % normalize
            parvals = tmpvals*(tmp(2)-tmp(1))+tmp(1);
        else
            error('Unidentified type of value distribution for param: "%s" given by %.2f.',Es.BFpar{ii},tmp(4));
        end;  
        
        if(ii==1)
            totvals=parvals;
        else
            if(abs(tmp(3))>0)
                tmpvals = repmat(totvals,curnum,1); % Make new points
                totvals = [tmpvals reshape(repmat(parvals,1,size(totvals,1))',curnum*size(totvals,1),1)];
            else
                len=size(totvals,1)/curnum;         % No new points
                totvals = [totvals reshape(repmat(parvals,1,len)',curnum*len,1)];
            end;
            
        end;
        
    end;
else % Or if we're in (1) format
    totvals = Es.BFrange;
end;

Es.BFparvals = totvals;
end
