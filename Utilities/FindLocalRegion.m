function [reg,relsz]=FindLocalRegion(regprm,Ps,Es,varargin)
% Find a local region with a given size ("volume") of the system
% reg=FindLocalRegion(regprm,Ps,Es)
% Where regprm = [sz x y], sz is the size of the system to choose,
% and (x,y) (or just x) is the location in the system
% all these can be "pixel/site" sized, or relative sized

% Update online if necessary
if(nargin>3) [~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:}); end;

if(regprm(1)<1)
    regprm(1)=ceil(regprm(1)*Ps.Nx*Ps.Ny);
end;

if(Ps.Ny>1) % for 2d system
    if(length(regprm)<3)
        if(regprm(2)<1)
            regprm(2)=ceil(regprm(2)*Ps.Nx*Ps.Ny);
        end;
        regprm(3)=ceil(regprm(2)/Ps.Nx);
        regprm(2)=mod(regprm(2)-1,Ps.Nx)+1;
    elseif(regprm(3)<1)
        regprm(3)=ceil(regprm(3)*Ps.Ny);
    end;
else % for 1d or network system
    if(regprm(2)<1)
        regprm(2)=ceil(regprm(2)*Ps.Nx);
    end;
end;

% assume periodic boundary conditions by default
if(~isfield(Ps,'Bc')) 
   Ps.Bc=0;
end;

if(isfield(Ps,'Net') && ~isempty(Ps.Net))
    inds=regprm(2); % a network structure
elseif(Ps.Ny==1)  % 1D case
    rad=regprm(1)/2-1;
    X  = meshgrid(-ceil(rad):ceil(rad),0);
    
    inds = X+regprm(2);
    if(Ps.Bc)
        inds((inds<1)|(inds>Ps.Nx))=[];
    else % periodic boundary conditions
        inds = mod(inds-1,Ps.Nx)+1;
    end;
else % 2D case
    rad=sqrt(regprm(1)/pi)-1;
    [X,Y]  = meshgrid(-ceil(rad):ceil(rad),-ceil(rad):ceil(rad));
    
    dist=sqrt(X.^2+Y.^2);
    X = X+regprm(2);
    Y = Y+regprm(3);
    
    if(Ps.Bc) 
        inside=~((X<1)|(X>Ps.Nx)|(Y<1)|(Y>Ps.Ny)|(dist>rad));
        inds=sub2ind([Ps.Nx Ps.Ny],X(inside),Y(inside));
    else % periodic boundary conditions
        X = mod(X-1,Ps.Nx)+1;
        Y = mod(Y-1,Ps.Ny)+1;
 
        inds=sub2ind([Ps.Nx Ps.Ny],X(:),Y(:));
        inds(dist>rad) = [];
    end;
end;

% Now, run recursively using nearest neighbor, to get exact number
nnsm=NeighborSM(1,Ps,Es);

% Create the region
reg=zeros(Ps.Nx*Ps.Ny,1);
reg(inds)=1;

counter = sum(reg);
maxiter = 1e4;
exitind = 1;
% Iteratively find neighboring sites
while counter<regprm(1) && exitind<maxiter  
    tmpreg = logical(nnsm*reg);
	tmpreg = tmpreg - reg.*tmpreg;
	reg = reg + tmpreg;
	counter = sum(reg ~= 0);
    exitind = exitind+1;
end;
if(exitind>=maxiter)
    warning('Reached maximum number of iterations while trying to find neighbors. Is network not connected?');
    while counter<regprm(1) && exitind<maxiter*2  % iteratively find neighboring sites
        loc=randi(Ps.Nx*Ps.Ny);
        if(~reg(loc))
            reg(loc)=1;
            counter=counter+1;
        end;
        exitind = exitind+1;
    end;
    if(exitind>=maxiter*2)
        warning('Could not find enough random locations either');
    end;
end;

if(counter>regprm(1))    % make sure we have exactly the right number of sites
	inds = find(tmpreg);
	reg(inds(1:(counter-regprm(1))))=0;     
end;

reg   = logical(reg);

relsz = regprm(1)/(Ps.Nx*Ps.Ny); % relative size of region found

end

