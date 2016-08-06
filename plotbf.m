function varargout=plotbf(bfs,Es,varargin)
% plot a bifurcation diagram
% plotbf(bfs,Es)
if(nargin<2)
	Es=struct();
end;

if(~isfield(Es,'ParmSpace'))
	Es.ParmSpace=0;
end;

if(~isfield(Es,'Balloon'))
	Es.Balloon=0;
end;

if(~isfield(Es,'BfColor'))
 	Es.BfColor=[0,0,0 ; hsv(9)];
end;

if(~isfield(Es,'BfPattern'))
 	Es.BfPattern=0;
end;

if(~isfield(Es,'BfFields'))
 	Es.BfFields=[1,2];
end;

if(~isfield(Es,'BfPhases'))
 	Es.BfPhases=[0,1];
end;

if(~isfield(Es,'BfStyle'))
	Es.BfStyle = ['-- ';'-  '];
end;

if(~isfield(Es,'Smooth'))
	Es.Smooth=0;
end;

mincount=3;
%ires=0.1;
small=0.00001;
ihold = ishold;

% Update online if necessary
if(nargin>3) [~,~,Es]=UpdateParameters([],[],Es,varargin{:}); end;

if(~iscell(bfs)) % put bfs into a cell array for convenience
	bfs={bfs};
end;

% Any negative BfFields value is considered to be counted from the end
Es.BfFields(Es.BfFields<0) = Es.BfFields(Es.BfFields<0)+size(bfs{1},2) +1;


%figure;
if(Es.ParmSpace&&(size(bfs,1)>1)&&(size(bfs,2)>1))
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
	DrawShape(temp,Es.BfFields,Es.Smooth,[Es.BfColor(1,:)]);
	ind=1;
	while(isempty(temp{ind}))
		ind=ind+1;
	end;
	
	for jj=1:length(Es.BfPhases)	
		hold on;
		DrawShape(temp,Es.BfFields,Es.Smooth,[Es.BfColor(jj+1,:) 0.5],[5 Es.BfPhases(jj)] );
	end;
	hold off;
	handle=gcf;
elseif(Es.Balloon & (length(bfs(:))>1))
	rind=1;
	Es.BfColor = [Es.BfColor  [1; 0.5*ones(size(Es.BfColor,1)-1,1)]];
	if(Es.Balloon==1)
		DrawShape(bfs,Es.BfFields,Es.Smooth,[Es.BfColor(1,:)]);
		ind=1;
		while(isempty(bfs{ind}))
			ind=ind+1;
		end;
		for jj=1:length(Es.BfPhases)
			hold on;
			DrawShape(bfs,Es.BfFields,Es.Smooth,Es.BfColor(jj+1,:),[size(bfs{ind},2) Es.BfPhases(jj)] );
		end;
		hold off;
		handle=gcf;
	else	% Draw balloon using bitmap image (and not the fill function)
		%Es.Balloon = [Es.Balloon(:) ; zeros(4,1)];
		im1 = flipdim(BalloonMat(bfs,Es.BfFields,Es.Balloon),1);
		if sum(Es.BfPattern(:))
			im1 = PatternBalloon(im1,Es.BfPattern(:,:,1));
		end;
		img = cat(3,im1*(1-Es.BfColor(1,1)),im1*(1-Es.BfColor(1,2)),im1*(1-Es.BfColor(1,3)))*Es.BfColor(1,4);
		ind = 1;
		while(isempty(bfs{ind}))
			ind=ind+1;
		end;
		for jj=1:length(Es.BfPhases)
			im2 = flipdim(BalloonMat(bfs,Es.BfFields,Es.Balloon,[size(bfs{ind},2) Es.BfPhases(jj)] ),1);
			imchange = repmat(im2>0,[1 1 3]);;
			if sum(Es.BfPattern(:))
				im2 = PatternBalloon(im2,Es.BfPattern(:,:,jj+1));
			end;
			
			imgcolor = cat(3,im2*(1-Es.BfColor(jj+1,1)),im2*(1-Es.BfColor(jj+1,2)),im2*(1-Es.BfColor(jj+1,3)));
			img(imchange) = img(imchange)*(1-Es.BfColor(jj+1,4)) + imgcolor(imchange)*Es.BfColor(jj+1,4);
			%DrawShape(bfs,Es.BfFields,Es.Smooth,[Es.BfColor(jj+1,:) 0.5],[size(bfs{ind},2) Es.BfPhases(jj)] );
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



function newbal = PatternBalloon(bal,pat)
repx = ceil(size(bal,1)/size(pat,1));
repy = ceil(size(bal,2)/size(pat,2));
tempbal = repmat(pat,repx,repy);
newbal = tempbal(1:size(bal,1),1:size(bal,2)).*bal;
end




