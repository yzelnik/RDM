function varargout=plotst(Vs,Ps,Es,stind,varargin)
% plot the state of some model
% plotst(Vs,Ps,Es,stind)

minst2colorlen = 10;
defst2colorlen = 64;
if(nargin<2)
	Ps=struct('Lx',size(Vs,1),'Nx',size(Vs,1),'Ny',1);
end;
if(nargin<3)
	Es=struct();
end;
if(nargin<4)
	stind=1;
end;

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if(~isfield(Es,'St1Color'))
	Es.St1Color = [jet(2) ;hsv(7)];
end;

if(~isfield(Es,'St2Color'))
	base = 1:-0.0138:0.13;
	Es.St2Color = [base; base/2+0.5; base]';
end;

if(~isfield(Es,'St2Sc'))
	Es.St2Sc = [];
end;

if(~isfield(Es,'St2Colorbar'))
	Es.St2Colorbar = 1;
end;

if(size(Es.St2Color,1)<minst2colorlen)
	Es.St2Color=interp1(0:(size(Es.St2Color,1)-1),Es.St2Color,(0:(defst2colorlen-1))/(defst2colorlen-1)*(size(Es.St2Color,1)-1),'cubic');
end;

% first check if this is a 1D plot
if((Ps.Nx==1) || (Ps.Ny==1))
	if(Ps.Nx==1)
		reallen=Ps.Ly;
	else	
		reallen=Ps.Lx;
	end;
	
	Ps.Nx=Ps.Nx*Ps.Ny;
	
    if (~isfield(Es,'Vind'))
	Es.Vind = 1:size(Vs,2);
    end;
    data=reshape(Vs(:,Es.Vind,stind),size(Vs,1),length(Es.Vind)*length(stind));
    set(0,'DefaultAxesColorOrder',[Es.St1Color(Es.Vind,:) ; Es.St1Color]);
    handle=plot((1:Ps.Nx)*(reallen/Ps.Nx),data);
    %if (isfield(Es,'Vind'))
%	set(0,'DefaultAxesColorOrder',[Es.St1Color(Es.Vind,:) ; Es.St1Color]);
%        plot((1:Ps.Nx*Ps.Ny)*(reallen/Ps.Nx*Ps.Ny),Vs(:,Es.Vind,stind));
%    else%
	%set(0,'DefaultAxesColorOrder',Es.St1Color);
        %plot((1:Ps.Nx*Ps.Ny)*(reallen/Ps.Nx*Ps.Ny),Vs(:,:,stind));
   % end
	xlim([0 Ps.Lx*Ps.Ny]);
else  % Assuming this is a 2D plot
	set(gcf,'Colormap',Es.St2Color);
	img = reshape(Vs(:,Es.Vind(1),stind),Ps.Nx,Ps.Ny)';
	if(isfield(Es,'St2Angle'))  % Rotate image if relevant
        img=imrotate(img,Es.St2Angle);
    end;
	if(~length(Es.St2Sc))	% Autoscale image?
		handle=imagesc(img);
	else
		handle=imagesc(img,Es.St2Sc);
	end;		
	set(gca,'XTickLabel',get(gca,'Xtick')*Ps.Lx/Ps.Nx); 
	set(gca,'YTickLabel',get(gca,'Ytick')*Ps.Ly/Ps.Ny);
	if(Es.St2Colorbar)      
		colorbar;
	end;
end

if(nargout>0)  % Only return a handle if one's requested.
    varargout{1}=handle;
end;


end
