function points=AnalyzeAutoStates(filename,Ps,Es,AnlFunc,ExtFunc,inds,updates,varargin)
% Read states from AUTO file and make an additional analysis on them
% points=AnalyzeAutoStates(filename,Ps,Es,AnlFunc,ExtFunc,inds,updates)
% Returns the points requested with the added results of the analysis
% AnlFunc is the function that makes the analysis, using ExtFunc if necessary 

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});

if(~isfield(Es,'AUTO'))
	Es.AUTO=[];
end;

% Read AUTO states file
[states,tags]=ReadAutoStates(sprintf('s.%s',filename),Ps,Es,0);

if strcmp(Es.AUTO,'notags')
	tmppoints=ReadAutoBif(sprintf('b.%s',filename),[0 inds(:)']);
	tags = find(tmppoints(:,1));
	tmppoints = tmppoints(:,2:end);
else
	% Read AUTO bifurcation file
	tmppoints=ReadAutoBif(sprintf('b.%s',filename),inds);
end;
%plotbf(tmppoints);
%plot(tmppoints(:,1),tmppoints(:,2))
%size(states)

if(~isfield(Es,'MultPeriod'))
	Es.MultPeriod=1;
end;

if(Es.MultPeriod(1)>1)
	%states = repmat(states,Es.MultPeriod(1),1);
	states=repmat(states,Es.MultPeriod(1),1);
    %basestate=states;
	%for ii=2:Es.MultPeriod
	%	states=[states; basestate];
	%end;
	%clear basestate;
	Ps.Nx=Ps.Nx*Es.MultPeriod(1);
end;
FirstLx=Ps.Lx;

points=[];
% Go through all states
for ii=1:length(tags)
	Ps.Lx=FirstLx;
	% Updating the fields from the bif file to the relevant state
	for jj=1:length(updates)   
		if(~isempty(updates{jj}))  % Is this field to be updated?
			if(ischar(updates{jj}))
				Ps.(updates{jj})=tmppoints(tags(ii),jj);
			else
				Ps.Ds(updates{jj})=tmppoints(tags(ii),jj);
			end;
		end;
	end;
	Ps.Lx=Ps.Lx*Es.MultPeriod;
	if(isfield(Ps,'SpaMat'))
	%	%disp('redoing mat')
		Ps.SpaMat = DervSM(2,Ps,Es);
	end;
	%ii
	%disp([ii Ps.Lx Ps.P Ps.Nx]);
	% Run the analysis function on the given state
    
    if isa(ExtFunc, 'function_handle')
		st=ExtFunc(states(:,:,ii),Ps,Es);
        %disp(sprintf('using %s!',func2str(ExtFunc)));
    else
        st=states(:,:,ii);
	end;
    res=AnlFunc(st,Ps,Es);
	% Add data to the returned array
    
    %plot(points(:,1),points(:,2));
    %pause;
	points = [points; tmppoints(tags(ii),:) res(:)'];
end;


