function Matt=DervSM(order,Ps,Es,varargin)
% Build a spatial matrix for a derivative of a given order
% Matt=DervSM(order,Ps,Es)

% Update online if necessary
if(nargin>3) [~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:}); end;

if(~isfield(Ps,'BC'))
   Ps.BC=0;
end;
len = Ps.Nx*Ps.Ny;	% Get size of the system (number of points)

% Is this a 1D system?
if((Ps.Nx==1) || (Ps.Ny==1))
	% Set correct resolution of the system
	if(Ps.Nx==1)
		res=Ps.Ly/Ps.Ny;
	else	
		res=Ps.Lx/Ps.Nx;
	end;

	switch order	% Pick relevant stencil for this order of derivative
   	  case 0
  		Matt=StencilToSM([0 1 0],Ps.Nx,Ps.Ny,Ps.BC);
  	  case 1
  		Matt=StencilToSM([1 0 -1],Ps.Nx,Ps.Ny,Ps.BC)/res/2;
  	  case 2
		Matt=StencilToSM([1 -2 1],Ps.Nx,Ps.Ny,Ps.BC)/res^2;
	  case 3
		Matt=StencilToSM([1 -2 0 2 -1],Ps.Nx,Ps.Ny,Ps.BC)/res^3/2;
	  case 4
		Matt=StencilToSM([1 -4 6 -4 1],Ps.Nx,Ps.Ny,Ps.BC)/res^4;
    	  otherwise
        	warning('High derivatives not supported (1D)');
		Matt=sparse(len,len);
	end;
else  % 2D system
	res=sqrt(Ps.Ly/Ps.Ny*Ps.Lx/Ps.Nx);
	switch order	% Pick relevant stencil for this order of derivative
   	  case 0
  		Matt=StencilToSM([0 0 0 ;0 1 0; 0 0 0],Ps.Nx,Ps.Ny,Ps.BC);
  	  case 1
  		Matt=StencilToSM([0 1 0;0 0 0;0 -1 0],Ps.Nx,Ps.Ny,Ps.BC)/res/2;
  	  case 2
		Matt=StencilToSM([0 1 0;1 -4 1; 0 1 0],Ps.Nx,Ps.Ny,Ps.BC)/res^2;
	  case 4
		Matt=StencilToSM([ 0 0 1 0 0; 0 2 -8 2 0; 1 -8 20 -8 1; 0 2 -8 2 0; 0 0 1 0 0],Ps.Nx,Ps.Ny,Ps.BC)/res^4;
    	  otherwise
        	warning('High derivatives not supported (2D)');
		Matt=sparse(len,len);
	end;
end;

