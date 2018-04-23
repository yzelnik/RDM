function Matt=DistSM(maxdist,Ps,Es,varargin)
% Build a distance spatial matrix (up to maxdist)
% Matt=DistSM(maxdist,Ps,Es)

% Update online if necessary
[~,Ps,~]=UpdateParameters([],Ps,Es,varargin{:});

if(~isfield(Ps,'Bc'))
   Ps.Bc=0;
end;

if((Ps.Nx==1) || (Ps.Ny==1))  
    [X,Y]  = meshgrid(-ceil(maxdist):ceil(maxdist),0);
else
    [X,Y]  = meshgrid(-ceil(maxdist):ceil(maxdist),-ceil(maxdist):ceil(maxdist));
end;

diststencil = sqrt(X.^2+Y.^2);
diststencil(diststencil>maxdist)=0;

Matt=StencilToSM(diststencil,Ps.Nx,Ps.Ny,-abs(Ps.Bc));  % make sure Bc is either 0 or -1


end

