function imm=BalloonMat(pnts,cols,dims,parse)
if(nargin<4)
	parse = 0;
end;

pnts=[pnts(:); {[]}];
imm = zeros(dims(1),dims(2));
wsmooth = 1;
ires = dims(3);

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
				newx=(min(xs):ires:max(xs))';
				bot=smooth(interp1(xs,mns,newx,'cubic'),ceil(wsmooth/ires));
				top=smooth(interp1(xs,mxs,newx,'cubic'),ceil(wsmooth/ires));	
				offset=round(newx(1)/dims(3));

				for ii=1:length(top)
%					[ii round(bot(ii)) round(top(ii))]
					imm(round(bot(ii)/dims(4)):round(top(ii)/dims(4)),ii+offset)=1;
				end;
				
		
			xs=[];
			mxs=[];
			mns=[];
		end;	
		ind=0;
	end;
end;

end
