function newbfs=ConcatBf(bfs,Es,varargin)
% Concat sets of bifurcation arrays, when possible
% newbfs=ConcatBf(bfs,Es)

distratio=15;   % How close should ends be to be considered concatable? (relative to avg dist)

if(nargin<2)
	Es=struct();
elseif nargin>2
    [~,~,Es]=UpdateParameters([],[],Es,varargin{:});
end;

if(~isfield(Es,'BfFields') && ~Es.BfFields)
 	Es.BfFields=[1,2];
end;

bfs = SortBf(bfs,Es);   % Make sure all bif arrays are properly sorted

if(~iscell(bfs))    % Wrap up in cell array form
	bfs={bfs};
end;

for ii=1:length(bfs)
    ends(ii*2+[-1 0],:)=[bfs{ii}(1,Es.BfFields);bfs{ii}(end,Es.BfFields)];
	avgdist(ii*2+[-1 0])=sqrt(mean(diag(pdist2(bfs{ii}(:,Es.BfFields),bfs{ii}(:,Es.BfFields)),1)));
end;
%ends
%avgdist


runflag=1;  % Run until there's no more good connections
while ((length(bfs)>1) && (runflag))
    %size(ends)
    endsdist=pdist2(ends,ends)./(avgdist'*avgdist);
    endsdist(logical(eye(length(avgdist))))=inf;
    endsdist(logical(diag([repmat([1 0],1,length(bfs)-1) 1],1)+diag([repmat([1 0],1,length(bfs)-1) 1],-1)))=inf;
    [minval,ind]=min(endsdist(:));
    %imagesc(endsdist)
    %ends
    %pause
    if(minval<distratio)
        
        % Find which 2 bif arrays connect "well"
        [a,b]=ind2sub(size(endsdist),ind);
        bfind=ceil([a b]/2);
        bfend=mod([a b],2);
        
        if(bfind(1)>bfind(2))   % sort these by their order
            bfind = bfind([2 1]);
            bfend = bfend([2 1]);
        end;
        
        % flip arrays if needed to allow concating 
        if(bfend(1)==1) 
            %disp('flip 1')
            bfs{bfind(1)}=flipdim(bfs{bfind(1)},1);
        end;
        if(bfend(2)==0)
            %disp('flip 2')
            bfs{bfind(2)}=flipdim(bfs{bfind(2)},1);
        end;
        %disp([minval bfind bfend])
        
        %pause
        % concat the two bif arrays into the first array
        bfs{bfind(1)}=[bfs{bfind(1)} ;bfs{bfind(2)}];
        % delete the other (old) bif array and its accompanying variables
        bfs(bfind(2))=[];
        ends(bfind(1)*2+[-1 0],:)=[bfs{bfind(1)}(1,Es.BfFields);bfs{bfind(1)}(end,Es.BfFields)];
        %ends(bfind(1)*2,:)=ends(bfind(2)*2,:);
        avgdist(bfind(1)*2+[-1 0])=sqrt(avgdist(bfind(1)*2)*avgdist(bfind(2)*2));
        ends(bfind(2)*2+[-1 0],:)=[];
        avgdist(bfind(2)*2+[-1 0])=[];
        
        %bfs(bfind(2):end-1)=bfs(bfind(2)+1:end);
        %ends(bfind(2)*2-1:end-2,:)=ends(bfind(2)*2+1:end,:);
        %avgdist(bfind(2)*2-1:end-2)=avgdist(bfind(2)*2+1:end);
    else
        runflag=0;
    end;
end;

newbfs=bfs;
%imagesc(endsdist)
end
