function newbf=RephaseBf(bfdst,bforg,Es,varargin)
% Add in "fake" phase values for one bif file (bfdst) based on another bif file (bforg)
% Using Es.fmod=1 can help results significantly (but may be slow)
% newbf=RephaseBf(bfdst,bforg,Es)

if(nargin<3)
	Es=struct();
elseif nargin>3 % Update online if necessary
    [~,~,Es]=UpdateParameters([],[],Es,varargin{:}); 
end;

if(~isfield(Es,'BfFields'))
 	Es.BfFields=[1,2];
end;

if(~isfield(Es,'fmod'))
   Es.fmod=0;
end;


cellflag=0;
if(~iscell(bfdst))    % Wrap up in cell array form
	bfdst={bfdst};
	cellflag=1;
end;
if(~iscell(bforg))    % Wrap up in cell array form
	bforg={bforg};
end;


if(length(bforg)~=length(bfdst))	% fix number of different bf arrays
	bforg = repmat(bforg(1),1,length(bfdst));
	warning('Writing-over/repeating bforg data');
end;

newbf={};
for ii=1:length(bforg)	% Go over each bf array
%	size(bforg{ii})
%size(bforg{ii}(:,Es.BfFields))

	org2dst_dist=pdist2(bforg{ii}(:,Es.BfFields),bfdst{ii}(:,Es.BfFields),'seuclidean');	% calculate dist between points
    %imagesc(org2dst_dist);
    %pause;
    if(Es.fmod==0)  % simple straight-forward method: find min distance per point
    	[minval,minind]=min(org2dst_dist);	% Find min dist per pair
        imagesc(org2dst_dist)
    else            % find consecutive (monotonic) points
        clear minind;
        dstind=1;
        lastval=inf;
        for orgind=1:(size(bforg{ii},1)-1)
            %disp([orgind dstind])
            %disp([org2dst_dist(orgind,dstind) lastval org2dst_dist(orgind+1,dstind)])
            while(org2dst_dist(orgind,dstind)<org2dst_dist(orgind+1,dstind) && dstind<(size(bfdst{ii},1)))
                minind(dstind)=orgind;
                %lastval=org2dst_dist(orgind,dstind);
                dstind=dstind+1;
            end;
            %dstind=dstind+1;
        end;
         minind(dstind:size(bfdst{ii},1))=orgind;
         %plot(minind);
        %pause;
       
    end;
    %size(minind)
    %[minval,minind]=min(org2dst_dist);	
    %plot(minind);
    %pause;
    %size(minind)

	newphs = bforg{ii}(minind,end);		% Get "fake" phases from min-pairs
	newbf{ii} = [bfdst{ii} newphs];		% Add fake phases to dstbf
end;

if(cellflag)		% unwrap bf, if it was recieved unwrapped
	newbf=newbf{1};
end;

end
