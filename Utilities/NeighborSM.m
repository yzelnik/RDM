function Matt=NeighborSM(order,Ps,Es,varargin)
% Build a nearest neighbors spatial matrix 
% Matt=NeighborSM(order,Ps,Es)

% Update online if necessary
[~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});

if(~isfield(Ps,'Bc'))
   Ps.Bc=0;
end;

if(order<1) 
    order = 1;
end;

if(isfield(Ps,'Net') && ~isempty(Ps.Net))
    Matt = Ps.Net;
else

    len = Ps.Nx*Ps.Ny;	% Get size of the system (number of points)

    if((Ps.Nx==1) || (Ps.Ny==1))    % 1D case
        if(order==1) 
            Matt=StencilToSM([1 0 1],Ps.Nx,Ps.Ny,Ps.Bc);                  
        else
            Matt=StencilToSM([1 1 0 1 1],Ps.Nx,Ps.Ny,Ps.Bc);                  
        end;
    else                            % 2D case
        if(order==1) 
            Matt=StencilToSM([0 1 0;1 0 1; 0 1 0],Ps.Nx,Ps.Ny,Ps.Bc);  
        else
            Matt=StencilToSM([1 1 1;1 0 1; 1 1 1],Ps.Nx,Ps.Ny,Ps.Bc);                  
        end;
    end;
end;
%Matt = logical(Matt);
end

