function varargout=plotbf(bfs,Es,varargin)
% plot a bifurcation diagram
% plotbf(bfs,Es)

if(nargin<2)
	Es=struct();
else
    if(mod(nargin,2)) error('No default extra-input exists for plotbf.'); end;
end;

% Update online if necessary
if(nargin>2) [~,~,Es]=UpdateParameters([],[],Es,varargin{:}); end;

% Put in some default values of Es
Es=InsertDefaultValues(Es,'PrmSpace',0,'BfBalloon',0,'BfColor',[0,0,0 ; hsv(9)],'BfPattern',0,'BfLineWidth',1);
Es=InsertDefaultValues(Es,'BfFields',[1,2],'BfPhases',[1,0],'BfStyle',['-  ';'-- '],'BfSmooth',0);

ihold = ishold; % is hold on/off activated?

if(~iscell(bfs)) % put bfs into a cell array for convenience
	bfs={bfs};
end;

ind=1;
    
% Any negative Es.BfFields value is considered to be counted from the end
bffields = Es.BfFields;
bffields(bffields<0) = bffields(bffields<0)+size(bfs{1},2) +1;

% Get color scheme
color=Es.BfColor(mod(ind-1,length(Es.BfColor))+1,1:3);
while(isempty(bfs{ind}))
	ind=ind+1;
end;
plot(bfs{ind}(1,bffields(1)),bfs{ind}(1,bffields(2)),'Color',color,'lineWidth',Es.BfLineWidth);
hold on;

% Run through all colors so legend will make sense with the colors used	
rind=1;
for ind=2:length(bfs(:))
	if(~isempty(bfs{ind}))
		color=Es.BfColor(mod(ind-1,length(Es.BfColor))+1,1:3);
		plot(bfs{ind}(1,bffields(1)),bfs{ind}(1,bffields(2)),'Color',color,'lineWidth',Es.BfLineWidth);
		rind=rind+1;
	end;	
end;

rind=1;
% And now actually plot	
for ind=1:length(bfs(:)) 
	if(~isempty(bfs{ind}))  % go over each cell that is not empty
        
    % Any negative Es.BfFields value is considered to be counted from the end
    bffields = Es.BfFields;
    bffields(bffields<0) = bffields(bffields<0)+size(bfs{ind},2) +1;

		% get data points and color to plot
        points=bfs{ind};
		color=Es.BfColor(mod(ind-1,length(Es.BfColor))+1,1:3);

        
        if(size(points,2)==1)   % allow for only a single column of data
			xx=1:size(points,1);	
			yy=points(:,1);
        else % Otherwise, take the two columns according to bffields
			xx=points(:,bffields(1));
			yy=points(:,bffields(2));
		end;
	
        if(Es.BfSmooth)
            ires=min(abs(diff(xx)))/10;
            newx=(min(xx):ires:max(xx))';
            newy=smooth(interp1(xx,yy,newx,'PCHIP'),Es.BfSmooth);
            xx=newx;
            yy=newy;
        end;
        
        % If more than two columns exist, assume that there is phase (i.e. stability) info
		if(size(points,2)>2)
            % Phase is given by Es.BfFields(3), or by default by the last column
            if(length(bffields)<3)
    			phase=points(:,end);
            else
                phase=points(:,bffields(3));
            end;
            if(Es.BfSmooth)
                phase=reshape(repmat(phase,1,10)',length(xx)+9,1);
            end;
			xs=xx(1);
			ys=yy(1);
            
            % go over each point in set of points
			for ii=2:size(xx,1)
                % concat with data so far
				xs=[xs; xx(ii)];
				ys=[ys; yy(ii)];
                % did we reach an end of a segment?
				if((phase(ii)~=phase(ii-1))||(ii==size(xx,1)))
					indphase=find(phase(ii-1)==Es.BfPhases);
					% figure out if this is a "known" phase
					if((isempty(indphase)) || (indphase(1)>size(Es.BfStyle,1)))
						plot(xs,ys,'Color',color,'lineWidth',Es.BfLineWidth); % unknown phase - simple plot
                    else
                        % we know the phase - plot accordingly
						plot(xs,ys,Es.BfStyle(indphase(1),:),'Color',color,'lineWidth',Es.BfLineWidth);
					end;
                    % start a new segment
					xs=xx(ii);
					ys=yy(ii);	
				
				end;
			end;
	
        else % no phase information exists, just simple plot
			handle=plot(xx,yy,Es.BfStyle(1,:),'Color',color,'lineWidth',Es.BfLineWidth);
		end;
	end;
	rind=rind+1;	
end;

if(ihold) % set back hold on/off to original setting
	hold on;
else
	hold off;
end;

if(nargout>0)  % Only return a handle if one's requested.
    varargout{1}=handle;
end;

end




