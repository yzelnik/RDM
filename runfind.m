function [StData,BfData]=runfind(Vs,Ps,Es,varargin)
% Repeatedly run Es.FindFunc (default is @runflow) 
% and get it close to Es.FindVal(1) (default is 0) 
% [StData,BfData]=runfind(Vs,Ps,Es)

% Default first extra input is for search-function to use
if(~mod(nargin,2)) varargin = ['Es.FindFunc' varargin]; end;
    
% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'FindFunc',@runflow,'FindVal',0,'FindOnlyBf',0);

Es.InitActive  = 0; % Allow states to be updated if necessary

if(length(Es.FindVal)<2)   % Define threshold value for search on the Bf-parameter(s)
    if(isfield(Es,'BfSmall') && ~isempty(Es.BfSmall))
        Es.FindVal(2)= Es.BfSmall;
    else
        Es.FindVal(2)= 1e-4;
    end;
end;
% Define threshold value for search for the test result itself
if((length(Es.FindVal)<3) || Es.FindVal(3)==0)
    Es.FindVal(3) = inf;
end;
if(Es.FindVal(2)==0)  
    Es.FindVal(2) = inf;
end;
if((Es.FindVal(2)==inf)&&(Es.FindVal(3)==inf))
	error('Both thresholds cannot be infinity');
end;
        
if(~iscell(Es.BfPrm))   % Wrap in cell array if not already in one
    Es.BfPrm={Es.BfPrm};
end;

% Load initial values from the available values at Ps
initvals = LoadParmList(Vs,Ps,Es); 

if(isnan(Es.FindVal(3)))
    % Run a heuristic minimization (no derivatives assumed)
    if(isfield(Es,'BfSmall') && Es.BfSmall>0)
        stepsize = Es.BfSmall*10;
    else
        stepsize = 1e-3;
    end;
    [finvals,testval,st] = HeuristicSearch(@(x) TestSearchFunc(Vs,Ps,Es,x,Es.FindVal(1),0),initvals,stepsize,Es.FindVal(2));
else
    % Define options for minimization
    opts = optimset('TolX',Es.FindVal(2),'TolFun',Es.FindVal(3));
    % Minimize test function
    finvals = fminsearch(@(x) TestSearchFunc(Vs,Ps,Es,x,Es.FindVal(1),1),initvals,opts);

    % run it again to get the state at the end location
    [~,testval,st]=TestSearchFunc(Vs,Ps,Es,finvals,Es.FindVal(1),1); 
end;

% Return results
BfData = [finvals testval];
if(Es.FindOnlyBf) % if we do not want to save the state-data
    StData = [];
else
    StData = st;
end;

end


%%%%%%%%%%%%%%%%%%%% AUX FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist,testval,st]=TestSearchFunc(Vs,Ps,Es,curpars,goodval,useabs) 

% Update parameters from curpars, and update everything else if needed
[Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,curpars);
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es);

% Run the system
[st,testval] = Es.FindFunc(Vs,Ps,Es);

%plotst(st,Ps,Es); title(testval); drawnow;
% Calculate the distance from the value we seek
if(useabs)
    dist = abs(goodval-testval(1));
else
    dist = goodval-testval(1);
end;

end
%%%%%%%% A heuristic search for non-continous values %%%%%%

function [finalpoint,testval,st]=HeuristicSearch(minfunc,initpoint,step,searchtol) 
if(length(initpoint(:))>1)
    error('Heuristic Search does not support more than one variable search!');
end;

initjumpnum=20;
jumps = reshape([1 -1]'*(2.^(0:initjumpnum-1)),initjumpnum*2,1)*step;

% First test
initval = minfunc(initpoint);
curval  = initval;

ind=1;
% Run over differnt jumps, look for a change in values
while((sign(initval)==sign(curval)) && ind<length(jumps))
    curpoint = initpoint + jumps(ind);
	curval   = minfunc(curpoint);
	ind=ind+1;
    %disp([ind curpoint curval])
end;

if(sign(initval)==sign(curval))
    warning('Search failed to find opposite values of objective function');
    finalpoint=NaN; 
    testval=[]; st=[];
else
	point1 = initpoint;
    point2 = curpoint;
    val1 = initval;
    val2 = curval;
    % Run until conditions are met
    while(abs(point1-point2)>searchtol)
        curpoint = (point1+point2)/2;  % binomial search
        curval   = minfunc(curpoint);  % run function at new position
        if(sign(curval)==sign(val1))
            val1   = curval;
            point1 = curpoint;
        elseif(sign(curval)==sign(val2))
            val2   = curval;
            point2 = curpoint;
        else  % we have 3 different signs of values (-,0,+)
            if(curval==0) || (val1==0)
                val2   = curval;
                point2 = curpoint;
            else
                val1   = curval;
                point1 = curpoint;
            end;    
        end;
        ind=ind+1;
        %disp([ind curpoint curval])
    end;

    finalpoint = curpoint;
    % Run at the result
    [~,testval,st] = minfunc(finalpoint);
end;

end

%%%%%%%%%%%  TO DELETE

%for jj=1:length(Es.BfPrm)
%	tmpval=curpars(jj);
%    if(isnumeric(Es.BfPrm{jj})) % Allow access to model-parameters by index
%        tmpfield=fieldnames(Ps);
%        Ps.(tmpfield{3+Es.BfPrm{jj}})=tmpval;
%    else
%        if(isempty(strfind(Es.BfPrm{jj},'.')) && ~strcmp(Es.BfPrm{jj},'Vs'))
%            Ps.(Es.BfPrm{jj}) = tmpval;  % for parameters in Ps (Prefereable)
%        else
%            eval(sprintf('%s=tmpval;',Es.BfPrm{jj}));
%        end;
%    end;
%end;

%disp([curpars dist])


%function initvals=GetInitVals(Vs,Ps,Es) 

%for jj=1:length(Es.BfPrm)
	
 %   if(isnumeric(Es.BfPrm{jj})) % Allow access to model-parameters by index
  %      tmpfield=fieldnames(Ps);
  %      tmpval=Ps.(tmpfield{3+Es.BfPrm{jj}});
  %  else
  %      if(isempty(strfind(Es.BfPrm{jj},'.')) && ~strcmp(Es.BfPrm{jj},'Vs'))
  %          tmpval=Ps.(Es.BfPrm{jj});  % for parameters in Ps (Prefereable)
  %      else
  %          eval(sprintf('tmpval=%s;',Es.BfPrm{jj}));
  %      end;
  %  end;
  %  initvals(jj)=tmpval;
%end;
   
%end


