function varargout=plotwf(Vs,Ps,Es,varargin)
% plot waterfall view (space-time plot)
% plotwf(Vs,Ps,Es)

if(nargin<2)
	Ps=struct('Lx',size(Vs,1),'Nx',size(Vs,1),'Ny',1);
end;
if(nargin<3)
	Es=struct();
else
    % Default first extra input is for the variable-indicator
    if(~mod(nargin,2)) varargin = ['Es.VarInd' varargin]; end;
end;

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);
% Put in some default values of Es
Es=InsertDefaultValues(Es,'VarInd',1,'SubPlot',[0 0.1 0.02 0.02 0.08 0.02]);

if(~isfield(Es,'St2Color') || ~Es.St2Color)
	base = 1:-0.0138:0.13;
	Es.St2Color = [base; base/2+0.5; base]';
end;


% first check if this is a 1D plot
if((Ps.Nx==1) || (Ps.Ny==1))
	if(Ps.Nx==1)
		reallen=Ps.Ly;
	else	
		reallen=Ps.Lx;
	end;
	
	Ps.Nx=Ps.Nx*Ps.Ny;
	set(gcf,'Colormap',Es.St2Color);	
	
    %length(Es.Frames)
    %size(Vs,3)
	if(isfield(Es,'Frames') && ~isempty(Es.Frames) && size(Vs,3)==length(Es.Frames) )
        yax = Es.Frames;
    else
        yax = 1:size(Vs,3);
    end;
    
    xax = (1:Ps.Nx)*(reallen/Ps.Nx);
    handle=imagesc(xax,yax,squeeze(Vs(:,Es.VarInd(1),:))');
    %set(gca,'XTickLabel',get(gca,'Xtick')*Ps.Lx/Ps.Nx); 
%	set(gca,'YTickLabel',get(gca,'Ytick')*Es.TimeDst/size(Vs,3));
	%plot((1:Ps.Nx*Ps.Ny)*(reallen/Ps.Nx*Ps.Ny),Vs(:,:,stind));
else  % Assuming this is a 2D plot
	clf
	frnum = size(Vs,3);
    if(Es.SubPlot(1))    % is there a pre-defined number of plots per row?
        subx = Es.SubPlot(1);
		suby = ceil(frnum/Es.SubPlot(1));
    else
    	if(frnum<6)
        	suby = 1;
        	subx = frnum;
    	elseif(frnum<9)
        	suby = 2;
        	subx = ceil(frnum/2);
        else
            subx = 6;
            suby = ceil(frnum/6);
        end;
    end;
	handle = tight_subplot(suby,subx,Es.SubPlot(2:3),Es.SubPlot(4:5),Es.SubPlot(6));
	for ii = 1:size(Vs,3); 
		axes(handle(ii)); 
		plotst(Vs(:,:,ii),Ps,Es);
	end;
	%warning('Waterfall view for 2D systems not implemented!');
	%imagesc(reshape(Vs(:,1),Ps.Nx,Ps.Ny),stind);
end

if(nargout>0)  % Only return a handle if one's requested.
    varargout{1}=handle;
end;

end
