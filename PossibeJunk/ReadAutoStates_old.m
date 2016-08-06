function [states indx]=ReadAutoStates(filename,Ps,Es,pnts,varargin)
% Read AUTO states file (usually s.name or fort.8)
% points=ReadAutoStates(filename,Ps,Es,inds)
% Returns a presentation of a set of states, identified by pnts

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});

% Initialization
if(nargin<4) 
	pnts=0;
end;
posthresh=-0.1;
abit = 0.001;

vrnum=Ps.Vnum;
fsize=Ps.Nx;
posflag=Es.posflag;

% Read the text File
data=dlmread(filename);

% rsz - number of sites, rlen - number of text lines in file per state, wid - number of columns (variables)
rsz=data(1,7);
rlen=data(1,9)+1;
wid=data(1,8);

% If points were not specified, return all of them
if (pnts==0)
	pnts = 1:round(length(data)/rlen);
end;


matt=[];
indx=[];
% go through all points requested, and read them
for pnt=pnts
	% Read this state
	tmpdata=(data((pnt-1)*rlen+2:(pnt-1)*rlen+rsz+1,1:wid));
	if(posflag)   % If we know values to be positive, then we can fix some AUTO bad writing of values
		pos=find((tmpdata<posthresh));
		tmpdata(pos)=10.^tmpdata(pos);
	end;
	% Add data to array
	matt=[matt tmpdata];
	indx=[indx data((pnt-1)*rlen+1,2)];
end;
%matt(end,:)=matt(end-1,:);
%matt(end,1:wid:end)=1.001;

% go through all points, and interpolate to a standard "grid" presentation
states=zeros(fsize,vrnum,length(pnts));
for ii=1:length(pnts)
	Xs=[0-abit ; matt(:,(ii-1)*wid+1) ; 1+abit ];
	Ys=[matt(1,(ii-1)*wid+(1:wid)) ; matt(:,(ii-1)*wid+(1:wid)) ; matt(end,(ii-1)*wid+(1:wid)) ];
	%[ size(matt(:,(ii-1)*wid+1)) size(matt(:,(ii-1)*wid+(1:wid))) size((0:(fsize-1))/fsize)]
	tmpmat=interp1(Xs,Ys,(0:(fsize-1))/fsize);
	states(:,:,ii)=tmpmat(:,1+(1:vrnum));
end;



