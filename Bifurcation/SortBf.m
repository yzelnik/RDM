function newbfs=SortBf(bfs,Es,varargin)
% sort a bifurcation array

%baseminnum = 100;

if(nargin<2)
	Es=struct();
elseif nargin>2
    [~,~,Es]=UpdateParameters([],[],Es,varargin{:});
end;

if(~isfield(Es,'BfFields'))
 	Es.BfFields=[1,2];
end;

if(~iscell(bfs))    % Wrap up in cell array form
	bfs={bfs};
end;

for ii=1:length(bfs)    % Go over the cells
    bf=bfs{ii};
    if(exist('pdist2','builtin')) % does function exist?
        dist=pdist2(bf(:,Es.BfFields),bf(:,Es.BfFields));
    else % or do it ourselves
        bfsz=size(bf,1);
        dist = zeros(bfsz);
        for ind=1:bfsz
            dist(:,ind) = sqrt(sum((bf(:,Es.BfFields)-repmat(bf(ind,Es.BfFields),bfsz,1)).^2,2));
        end;
    end;
    
    %imagesc(dist)
    dist(logical(eye(size(bf,1))))=inf;
    [~,ind]=min(dist(:));
    [a,b]=ind2sub(size(dist),ind);
    
    chain=[a b];
    for jj=1:size(bf,1)-2
        dist(chain(1),chain(end))=inf;
        dist(chain(end),chain(1))=inf;
        [two,ind]=min(dist(:,[chain(1) chain(end)]));
        if(two(1)<two(2))
            %([jj 1 ind(1,1)])
            %ind
            dist(chain(1),:)=inf;
            chain = [ind(1,1) chain];
            
        else
            %([jj 2 ind(1,2)])
            dist(chain(end),:)=inf;
            chain = [chain ind(1,2)];
        end;
        %chain
    end;
    newbfs{ii}=bfs{ii}(chain,:);
    %basenum=max(baseminnum,size(bf,1)*4);
    %bases=interp1(0:size(bf,1)-1,bf(:,Es.BfFields),(1:basenum)*((size(bf,1)-1)/basenum));
    %difmat = squeeze(sum((repmat(bf(:,Es.BfFields),[1 1 basenum])-shiftdim(repmat(bases',[1 1 size(bf,1)]),2)).^2,2));
    %[~,minind]=min(difmat,[],2);
    %[~,newind] = sort (minind);
    %newbfs{ii} = bf(newind,:);
    %plot(newind)
   % size(repmat(bf(:,Es.BfFields),[1 1 basenum]))
   % size(shiftdim(repmat(bases',[1 1 size(bf,1)]),2))
    %size(sum(repmat(bf(:,Es.BfFields),[1 1 basenum])-shiftdim(repmat(bases',[1 1 size(bf,1)]))))
    %imagesc(log10(difmat));
    %imagesc(bf(:,Es.BfFields)*bases');
    %plot(bases);
    %for jj=1:length(Es.BfFields)
        
        %bases(:,jj)=interp1(1:size(bf,1),bf(:,Es.BfFields(jj)),(1:basenum)*(size(bf,1)/basenum));
    %end;
end;


if(length(newbfs)==1)   % Unwarp in necessary
    newbfs=newbfs{1};
end;

end

