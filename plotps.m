function varargout=plotps(bfs,Es,varargin)
% plot a parameter-space
% plotps(bfs,Es)
% This assumes Es.BfFields has 3 values,
% for the 2 axes and where the values to show are

griderr = 1e-10; % how much error can a grid take?

if(nargin<2)
	Es=struct();
else
    if(mod(nargin,2)) error('No default extra-input exists for plotps.'); end;
end;

% Update online if necessary
if(nargin>2) [~,~,Es]=UpdateParameters([],[],Es,varargin{:}); end;

% Put in some default values of Es
Es=InsertDefaultValues(Es,'BfFields',[1,2,3],'BfPhases',[0,1],'BfFilter',[],'BfSmooth',0);

if(isempty(bfs)) % make sure there's data
    error('No data to plot.');
end;

% filter according to Es.BfFilter if relevant
if(~isempty(Es.BfFilter))
    if(size(Es.BfFilter,2)<2) % set to horizontal vector if needbe
        Es.BfFilter=Es.BfFilter';
    end;
    if(size(Es.BfFilter,2)<3) % padd in zeros if need-be
        Es.BfFilter(1,3)=0;
    end;
    for ii=1:size(Es.BfFilter,1) % go over each field to filter
        tmp = abs(bfs(:,Es.BfFilter(ii,1))-Es.BfFilter(ii,2))<=Es.BfFilter(ii,3);
        bfs=bfs(tmp,:);
    end;
end;

% make sure there's enough data and fields
if(isempty(bfs))
    error('No data to plot after filtering.');
end;
if(length(Es.BfFields)<3)
    error('3 fields needed to make a parameter-space plot');
end;

% figure out the 2 possible axes, and if this is a grid
yax=unique(bfs(:,Es.BfFields(1)));
xax=unique(bfs(:,Es.BfFields(2)));
isgrid = griderr>(sum(abs(diff(diff(xax))))+sum(abs(diff(diff(yax)))));

% if this is a grid, plot it out
if(isgrid)
    if(size(bfs,1)==(length(xax)*length(yax)))
        grid = reshape(bfs,length(yax),length(xax),size(bfs,2));
    else
        % if things do not match up exactly, go pixel by pixel
        grid = zeros(length(yax),length(xax),size(bfs,2));
        for ii=1:length(yax) % go over axis 1
            tmp1 = bfs(bfs(:,Es.BfFields(1))==yax(ii),:);
            for jj=1:length(xax) % go over axis 2
                tmp2=tmp1(tmp1(:,Es.BfFields(2))==xax(jj),:);
                grid(ii,jj,:)=mean(tmp2,1); % calculate average
            end;
        end;
    end;
    % plot it out
    handle=imagesc(yax,xax,grid(:,:,Es.BfFields(3))');
    axis xy;
else % just make a scatter plot
    handle=scatter(bfs(:,Es.BfFields(1)),bfs(:,Es.BfFields(2)),20,bfs(:,Es.BfFields(3)));
end;
%axis([min(bfs(:,Es.BfFields(1))) max(bfs(:,Es.BfFields(1))) min(bfs(:,Es.BfFields(2))) max(bfs(:,Es.BfFields(2)))]);

if(nargout>0)  % Only return a handle if one's requested.
    varargout{1}=handle;
end;

end





% this is some old code for figuring out a grid...
%[sum(abs(diff(diff(xax)))) sum(abs(diff(diff(yax))))]

% Find the first point where the parameter Es.BfFields(1) changes direction
%tmp=diff(bfs(:,Es.BfFields(1)));
%switchpoints = find(~(sign(tmp(1))==sign(tmp)));
% Using this, try to see what is the size of this as a grid, if it is one
%leny = switchpoints(1);
%lenx = size(bfs,1)/leny;

%if(mod(lenx,1))
%    lenx=0;  % if the numbers don't add up, this is not a grid
%else
    % if the number of entries are correct, check if indeed this is a grid
%    grid = reshape(bfs,leny,lenx,size(bfs,2));
    
    % Check if the whole grid is completely consistent
 %   griderr = sum(sum(abs(diff(grid(:,:,Es.BfFields(1))'))))+sum(sum(abs(diff(grid(:,:,Es.BfFields(2))))));
 %   if(griderr)
 %       leny=0; % not a good grid
 %   else
        % if things are consistent, get the axis correctly
  %      yax  = bfs(1:leny,Es.BfFields(1));
  %      xax  = bfs(1:leny:end,Es.BfFields(2));
  %  end;
    
%end;



