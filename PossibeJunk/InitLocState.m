function VsOut=InitLocState(Vs,Ps,Es,repnum,phases,uniform,finsize,varargin)
% Initiate a localized state
% VsOut=InitLocState(Vs,Ps,Es,repnum,phases,uniform,finsize)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% First check if this is a 1D system
if((Ps.Nx==1) || (Ps.Ny==1))
	Ps.Nx=Ps.Nx*Ps.Ny;
		
	% Fill out the "empty" phases with zeros of none we given
	phases = [phases(:);0;0;0];
	
	% Change the phase of the initial (periodic) solution used
	locst = [Vs(round(phases(1)*Ps.Nx)+1:end,:,1) ; Vs(1:round(phases(1)*Ps.Nx),:,1)];
	
	% Repeat the periodic solution repnum times
	finst=[]; 	
	for ii=1:repnum 
		finst=[finst ; locst]; 
	end;
	
	% Fill out the solution with uniform values
	uniflen = finsize-size(finst,1);
	finst = [ones(floor(uniflen/2),1)*uniform(:)' ; finst ; ones(ceil(uniflen/2),1)*uniform(:)'];
	
	% Change the overall phase of the final solution	
	VsOut = [finst(round(phases(2)*finsize)+1:end,:,1) ; finst(1:round(phases(2)*finsize),:,1)];

else  % Assuming this is a 2D system
	% NOTE! the 2D version "uses" the input paramaters with different meaning. 
	% This should be changed.
	Vs = reshape(Vs,Ps.Nx,Ps.Ny,Ps.Vnum);
	VsOut = ones(finsize(1)*finsize(2),1)*uniform(:)';
	VsOut = reshape(VsOut,finsize(1),finsize(2),size(VsOut,2));

	iy = phases(1); 
	ix = phases(2); 
	rad = repnum; 
	for ii=1:finsize(1) 
		for jj=1:finsize(2)
			if(sqrt((ix-ii)^2+(iy-jj)^2)<rad) 
				VsOut(ii,jj,:)=Vs(ii,jj,:);
			end;  
		end; 
	end;
	VsOut = reshape(VsOut,finsize(1)*finsize(2),size(VsOut,3));
end

end
