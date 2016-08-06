function [num,stats]=T_CountRegions(Vs,Ps,Es,varargin)
% Counts the number of different regions in the given state
% stat gives, per region type: [sum sumvar sz szvar]
% where sum and sz are the average (per region) of the total sum of the region and the size of it respectively
% and sumvar and szvar are the respective variances of these values 

% data is assumed to be in the form (x*y*z) where x&y are spatial and z is the different state variables
% varparms = [vnum vindex dim1 dim2] 
% vnum is the number of state variables and vindex chooses which one to do the statistics on
% dim1 and dim2 give the x/y size of the data
% thresh gives the minimum value of the relevant variable to still count within the region

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

dim1 = Ps.Nx;
dim2 = Ps.Ny;
vnum = Ps.Vnum;
vindex = Es.Vind(1);

vmax = max(Vs(:,vindex,1));
vmin = min(Vs(:,vindex,1));

if((vmax-vmin)<Es.STsmall)
	num =[1 0];
	stats=zeros(1,4);
else

thresh = (vmax-vmin)/2;

for ind=1:2
	if(ind==1)
		sdata = reshape(Vs(:,vindex,1),dim1,dim2)-vmin;
	else
		sdata = vmax-reshape(Vs(:,vindex,1),dim1,dim2);
	end;
	%datas=zeros(1000,length(sdata(:)));
	regs = zeros(1000,2);
	regnum = 0;
	[aa,bb] = max(sdata(:));
	while(aa>thresh)
		[xx,yy] = ind2sub(size(sdata),bb);
		regnum=regnum+1;
		sdata(xx,yy)=-1;
		regs(regnum,1)=regs(regnum,1)+aa;	
		regs(regnum,2)=regs(regnum,2)+1;
		flag=1;
		while(flag)
			flag=0;
			for ii=1:dim1 for jj=1:dim2
				if(sdata(ii,jj)==-1)
					sdata(ii,jj)=-2;					
					for x=ii-1:ii+1 for y=jj-1:jj+1
						if(x<1) rx=dim1; else if(x>dim1) rx=1; else rx=x; end; end;
						if(y<1) ry=dim2; else if(y>dim2) ry=1; else ry=y; end; end;
						
						if(sdata(rx,ry)>thresh)
							regs(regnum,1)=regs(regnum,1)+sdata(rx,ry);
							regs(regnum,2)=regs(regnum,2)+1;
							sdata(rx,ry)=-1;
							flag=1;
						end;
				end; end;
				end;
			end; end;
		end;
		[aa,bb] = max(sdata(:));
		%datas(regnum,:)=sdata(:);
	end;
	regs=regs(1:regnum,:);
    
	stats(ind*4-3:ind*4) = [mean(regs(:,1)) sqrt(mean((regs(:,1)-mean(regs(:,1))).^2)) mean(regs(:,2)) sqrt(mean((regs(:,2)-mean(regs(:,2))).^2)) ];
	num(ind) = regnum;
end;

end;
%results=regs(1:regnum,:);
%datas=datas(1:regnum,:);

end
