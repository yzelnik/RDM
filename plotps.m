function varargout=plotps(bfs,Es,varargin)
% plot a parameter-space
% plotps(bfs,Es)

if(nargin<2)
	Es=struct();
else
    if(mod(nargin,2)) error('No default extra-input exists for plotps.'); end;
end;

% Update online if necessary
if(nargin>2) [~,~,Es]=UpdateParameters([],[],Es,varargin{:}); end;

% Put in some default values of Es
%Es=InsertDefaultValues(Es,'PrmSpace',0,'BfBalloon',0,'BfColor',[0,0,0 ; hsv(9)],'BfPattern',0);
Es=InsertDefaultValues(Es,'BfFields',[1,2,3],'BfPhases',[0,1],'BfFilter',[],'BfSmooth',0);

if(isempty(bfs))
    error('No data to plot.');
end;

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

if(isempty(bfs))
    error('No data to plot after filtering.');
end;

if(length(Es.BfFields)<3)
    error('3 fields needed to make a parameter-space plot');
end;

% Find the first point where the parameter Es.BfFields(1) changes direction
tmp=diff(bfs(:,Es.BfFields(1)));
switchpoints = find(~(sign(tmp(1))==sign(tmp)));
% Using this, try to see what is the size of this as a grid, if it is one
leny = switchpoints(1);
lenx = size(bfs,1)/leny;
if(mod(lenx,1))
    lenx=0;  % if the numbers don't add up, this is not a grid
else
    % if the number of entries are correct, check if indeed this is a grid
    grid = reshape(bfs,leny,lenx,size(bfs,2));
    
    % Check if the whole grid is completely consistent
    griderr = sum(sum(abs(diff(grid(:,:,Es.BfFields(1))'))))+sum(sum(abs(diff(grid(:,:,Es.BfFields(2))))));
    if(griderr)
        leny=0; % not a good grid
    else
        % if things are consistent, get the axis correctly
        yax  = bfs(1:leny,Es.BfFields(1));
        xax  = bfs(1:leny:end,Es.BfFields(2));
    end;
    
end;

if(lenx*leny)
    handle=imagesc(yax,xax,grid(:,:,Es.BfFields(3))');
    axis xy;
else
    handle=scatter(bfs(:,Es.BfFields(1)),bfs(:,Es.BfFields(2)),20,bfs(:,Es.BfFields(3)));
end;
%axis([min(bfs(:,Es.BfFields(1))) max(bfs(:,Es.BfFields(1))) min(bfs(:,Es.BfFields(2))) max(bfs(:,Es.BfFields(2)))]);

if(nargout>0)  % Only return a handle if one's requested.
    varargout{1}=handle;
end;

end




function dummydunc()

mincount=3;
small=0.00001;
ihold = ishold;

if(~iscell(bfs)) % put bfs into a cell array for convenience
	bfs={bfs};
end;

% Any negative BfFields value is considered to be counted from the end
Es.BfFields(Es.BfFields<0) = Es.BfFields(Es.BfFields<0)+size(bfs{1},2) +1;
size(bfs)
%figure;
if(Es.PrmSpace&&(size(bfs,1)>1)&&(size(bfs,2)>1))
    disp(111)
	temp={};
	for ii=1:size(bfs,1)
		temp{ii}=[];
		for jj=1:size(bfs,2)
			if(~isempty(bfs{ii,jj}))
				unq=unique(bfs{ii,jj}(:,end));
				newvals=[];	
				for kk=1:length(unq)
					if(length(nonzeros(bfs{ii,jj}(:,end)==unq(kk)))>mincount)
						newv = [bfs{ii,jj}(1,1:end-1) unq(kk)];
						
						newv(Es.BfFields(2)) = newv(Es.BfFields(2)) + unq(kk)*small;
						newvals=[newvals; newv];
					end;			
				end;
				%temp{ii}=[temp{ii}; [ones(size(unq))*bfs{ii,jj}(1,Es.BfFields(1)) unq*small+bfs{ii,jj}(1,Es.BfFields(2)) unq]];	
				temp{ii} = [temp{ii}; newvals];	
			end;
		end;
	end;
%	temp
	DrawShape(temp,Es.BfFields,Es.BfSmooth,[Es.BfColor(1,:)]);
	ind=1;
	while(isempty(temp{ind}))
		ind=ind+1;
	end;
	
	for jj=1:length(Es.BfPhases)	
		hold on;
		DrawShape(temp,Es.BfFields,Es.BfSmooth,[Es.BfColor(jj+1,:) 0.5],[5 Es.BfPhases(jj)] );
	end;
	hold off;
	handle=gcf;
elseif(Es.BfBalloon & (length(bfs(:))>1))
	rind=1;
	Es.BfColor = [Es.BfColor  [1; 0.5*ones(size(Es.BfColor,1)-1,1)]];
	if(Es.BfBalloon==1)
		DrawShape(bfs,Es.BfFields,Es.BfSmooth,[Es.BfColor(1,:)]);
		ind=1;
		while(isempty(bfs{ind}))
			ind=ind+1;
		end;
		for jj=1:length(Es.BfPhases)
			hold on;
			DrawShape(bfs,Es.BfFields,Es.BfSmooth,Es.BfColor(jj+1,:),[size(bfs{ind},2) Es.BfPhases(jj)] );
		end;
		hold off;
		handle=gcf;
	else	% Draw balloon using bitmap image (and not the fill function)
		%Es.BfBalloon = [Es.BfBalloon(:) ; zeros(4,1)];
		im1 = flipdim(BalloonMat(bfs,Es.BfFields,Es.BfBalloon),1);
		if sum(Es.BfPattern(:))
			im1 = PatternBalloon(im1,Es.BfPattern(:,:,1));
		end;
		img = cat(3,im1*(1-Es.BfColor(1,1)),im1*(1-Es.BfColor(1,2)),im1*(1-Es.BfColor(1,3)))*Es.BfColor(1,4);
		ind = 1;
		while(isempty(bfs{ind}))
			ind=ind+1;
		end;
		for jj=1:length(Es.BfPhases)
			im2 = flipdim(BalloonMat(bfs,Es.BfFields,Es.BfBalloon,[size(bfs{ind},2) Es.BfPhases(jj)] ),1);
			imchange = repmat(im2>0,[1 1 3]);;
			if sum(Es.BfPattern(:))
				im2 = PatternBalloon(im2,Es.BfPattern(:,:,jj+1));
			end;
			
			imgcolor = cat(3,im2*(1-Es.BfColor(jj+1,1)),im2*(1-Es.BfColor(jj+1,2)),im2*(1-Es.BfColor(jj+1,3)));
			img(imchange) = img(imchange)*(1-Es.BfColor(jj+1,4)) + imgcolor(imchange)*Es.BfColor(jj+1,4);
			%DrawShape(bfs,Es.BfFields,Es.BfSmooth,[Es.BfColor(jj+1,:) 0.5],[size(bfs{ind},2) Es.BfPhases(jj)] );
		end;
		img = 1-img;
		handle=image(img);
		%image(img);
		%xtl=get(gca,'XTickLabel');
		%set(gca,'XTickLabel',xtl);
	
	end;

else
	ind=1;
    
	color=Es.BfColor(mod(ind-1,length(Es.BfColor))+1,1:3);
	while(isempty(bfs{ind}))
		ind=ind+1;
	end;
	plot(bfs{ind}(1,Es.BfFields(1)),bfs{ind}(1,Es.BfFields(2)),'Color',color);
	hold on

	% Run through all colors so legend will make sense with the colors used	
	rind=1;
	for ind=2:length(bfs(:))
		if(~isempty(bfs{ind}))
			color=Es.BfColor(mod(ind-1,length(Es.BfColor))+1,1:3);
			plot(bfs{ind}(1,Es.BfFields(1)),bfs{ind}(1,Es.BfFields(2)),'Color',color);
			rind=rind+1;
		end;	
	end;
	rind=1;
	% And now actually plot	
	for ind=1:length(bfs(:))
		if(~isempty(bfs{ind}))
			points=bfs{ind};
			color=Es.BfColor(mod(ind-1,length(Es.BfColor))+1,1:3);

			if(size(points,2)==1)
				xx=1:size(points,1);	
				yy=points(:,1);
			else
				xx=points(:,Es.BfFields(1));
				yy=points(:,Es.BfFields(2));
			end;
		
			if(size(points,2)>2)
				phase=points(:,end);
				xs=xx(1);
				ys=yy(1);
				
				for ii=2:size(points,1)
					xs=[xs; xx(ii)];
					ys=[ys; yy(ii)];
					if((phase(ii)~=phase(ii-1))||(ii==size(points,1)))
						indphase=find(phase(ii-1)==Es.BfPhases);
						%disp(indphase);
						%disp(ii)
						
						if((isempty(indphase)) | (indphase(1)>size(Es.BfStyle,1)))
							plot(xs,ys,'Color',color);
						else
							plot(xs,ys,Es.BfStyle(indphase(1),:),'Color',color);
						%	Es.BfStyle(indphase(1),:)
						%	size(xs)
						end;
						xs=xx(ii);
						ys=yy(ii);	
					
					end;
				end;
		
			else
				handle=plot(xx,yy,Es.BfStyle(1,:),'Color',color);
			end;
		end;
		rind=rind+1;	
	end;
	hold off
end;

if(ihold)
	hold on;
else
	hold off;
end;

if(nargout>0)  % Only return a handle if one's requested.
    varargout{1}=handle;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newbal = PatternBalloon(bal,pat)
repx = ceil(size(bal,1)/size(pat,1));
repy = ceil(size(bal,2)/size(pat,2));
tempbal = repmat(pat,repx,repy);
newbal = tempbal(1:size(bal,1),1:size(bal,2)).*bal;
end




