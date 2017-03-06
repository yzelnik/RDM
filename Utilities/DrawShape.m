function DrawShape(pnts,cols,wsmooth,color,parse,xy)
% Draw a shape
% par = scanned parameter; U = result (variable at steady state); K = wavenumber;
if(nargin<5)
	parse = 0;
end;
if(nargin<6)
	xy = 0;
end;

ires=10;
pnts=[pnts(:); {[]}];
color=[color(:) ;1]';

ind=0;
for ii=1:length(pnts)
	temp=pnts{ii};
	if((length(parse)>1) & (size(temp)>0))
			temp=temp(temp(:,parse(1))==parse(2),:);
	end;
	if(size(temp)>0)
		ind=ind+1;
		
		xs(ind)=temp(1,cols(1)); 			
		mxs(ind)=max(temp(:,cols(2)));
		mns(ind)=min(temp(:,cols(2))); 
	elseif (ind>0)

		if(ind>1)
			if(wsmooth(1)>0)
				newx=(min(xs):ires:max(xs))';
				bot=smooth(interp1(xs,mns,newx,'PCHIP'),ceil(wsmooth/ires));
				top=smooth(interp1(xs,mxs,newx,'PCHIP'),ceil(wsmooth/ires));	
            elseif(wsmooth(1)<0)
                newx=xs';
	            bot=cleancurve(mns',-wsmooth(1));
                top=cleancurve(mxs',-wsmooth(1));
                if(length(wsmooth)>1)
                    bot=cleancurve(bot,-wsmooth(2));
                    top=cleancurve(top,-wsmooth(2));
                end;
	        else
				newx=xs';
				bot=mns';
				top=mxs';
							end;
			%size(newx)
			%size(bot)
						%disp([newx(1) newx(end); top(1) top(end); bot(1) bot(end)]);
			%[newx top bot]
			if(xy)
                fill([top ; bot(end:-1:1)],[newx ; newx(end:-1:1)],color(1:3), 'FaceAlpha', color(4),'EdgeColor','none');			
            else
                fill([newx ; newx(end:-1:1)],[top ; bot(end:-1:1)],color(1:3), 'FaceAlpha', color(4),'EdgeColor','none');			
            end;
            hold on;
			xs=[];
			mxs=[];
			mns=[];
		end;	
		ind=0;
	end;
end;
hold off;
end


function curve=cleancurve(curve,jump)

	tmplist  = find(abs(diff(curve(:)))>jump);
	finlist  = tmplist(diff(tmplist)==1);
    jumpflag = abs(curve(finlist)-curve(finlist+2))<min([abs(curve(finlist)-curve(finlist+1)) abs(curve(finlist+1)-curve(finlist+2))],[],2);
    curve(finlist(jumpflag)+1)=(curve(finlist(jumpflag))+curve(finlist(jumpflag)+2))/2;         

end