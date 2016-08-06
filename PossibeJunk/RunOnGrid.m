function results=RunOnGrid(Vs,Ps,Es,IntFunc,AnlFunc,pname,pgrid,varargin)
% Run an analyis function on many different points (in parameter space) in some grid
% results=RunOnGrid(Vs,Ps,Es,IntFunc,AnlFunc,pname,pgrid,varargin)
% pname is a cell array of parameter names, pgrid is a cell array of a vector of points of values for this parameter
% AnlFunc is the function that makes the analysis, more parameters should follow within varargin

% Note - NO online update

% Number of different parameters in grid
plen=length(pname);

% The different indices
inds=ones(plen,1);

% current parameter values (for writing the outpout)
modeinfo=inds;

% Find the max value each index gets
maxinds=inds;
for ii=1:plen
	maxinds(ii)=length(pgrid{ii}(:));
end;
% The total values of the indices (used in the while loop)
totlen=sum(maxinds);

results=[];
inds(1)=0;
% The strange way of making several for loops together
while(sum(inds)<totlen)
	% Push indices forward
	inds(1)=inds(1)+1;
	
	% If an index reaches the end, put it back to the start, and push the next one a bit
	for ii=1:plen
		if(inds(ii)>maxinds(ii))
			inds(ii)=1;
			inds(ii+1)=inds(ii+1)+1;
		end;
	end;
	
	% Update parameters, both in the "Ps"- for the actual run, and in the "modeinfo" for the output
	for ii=1:plen   
		Ps.(pname{ii})=pgrid{ii}(inds(ii));
		modeinfo(ii)=pgrid{ii}(inds(ii));
	end;	
	out=RunToSS(Vs,Ps,Es,IntFunc);
	% Run the function
	tmp=AnlFunc(out,Ps,Es,varargin{:});
	% Summarize the results
	results = [results; modeinfo' tmp(:)'];
end;


