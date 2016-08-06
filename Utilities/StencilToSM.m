function Matt=StencilToSM(stencil,nx,ny,BC)
% Use a stencil to build a spatial matrix
% Matt=StencilToSM(stencil,nx,ny,BC)
% nx,ny give the number of sites in the system
% The stencil is assumed to be smaller than (nx,ny) with odd number of sites. 
% BC gives the Boundary Conditions, 0=Periodic, 1=Neumann, -1=Dirichlet
% If BC is not specifed, default = 0 (periodic)

if(nargin<4)
    BC=0;   % by default boundary conditions (BC) are periodic
end;
if(all(~BC))   
    BCflag=0;   % All BC's are periodic
else
    BCflag=1;   % some/all BC are NOT periodic
end;

if ~prod(mod(size(stencil),2))
    error('Illegal stencil: size must be odd in each axis.');
end;

% Initialize
Matt=sparse(nx*ny,nx*ny);

% if this is 1D, make sure things are set up correctly
if((nx==1) || (ny==1))
	stencil=stencil(:);
	nx=nx*ny;
	ny=1;
end;

szx = (size(stencil,1)-1)/2;
szy = (size(stencil,2)-1)/2;
% Go over all the stencil 
for ix=-szx:szx
	for iy=-szy:szy
		if(stencil(ix+szx+1,iy+szy+1))
			% Find the absolute position of a stencil site, considering the full size of the system (nx,ny)
			sx = mod(ix,nx);
			sy = mod(iy,ny);
			% Build a vector containing the values along the diagonal
			base = repmat([ones(nx-sx,1) ;zeros(sx,1)],ny,1);

            % Build a matrix with the diagonal of base (cyclic, so that there might be 2 diagonals as a result)
            tmpmat = spdiags(repmat(base,1,2),[ -(sx+sy*nx) nx*ny-(sx+sy*nx)],nx*ny,nx*ny);
            % more efficient than:   tmpmat = sparse(circshift(diag(base),sx+sy*nx));
			% If the lesser (small) (X) axis is not zero, than an additional diagonal is needed (from the cyclic in the small axis)  
			if(sx~=0)
				% Similarly to before, add another diagonal (or effectively 2) is the same way
                tmpmat = tmpmat + spdiags(repmat(1-base,1,2),[ -(sx+sy*nx-nx) nx*ny-(sx+sy*nx-nx)],nx*ny,nx*ny);
                % more efficient than:   tmpmat = tmpmat + sparse(circshift(diag(1-base),sx+sy*nx-nx));	
            end;
            if(BCflag && (ix||iy))
                changevec = CorrectBC(ix,iy,nx,ny,BC);
                for ind=1:size(changevec,3)   % go over all changes (first one is deleting the original periodic sites)
                    BCchangemat = sparse(changevec(:,1,ind),changevec(:,2,ind),changevec(:,3,ind),nx*ny,nx*ny);
                    tmpmat = tmpmat+BCchangemat;
                end;
                %zeros(size(tmpmat)); 
                %BCchangemat(sub2ind([nx*ny,nx*ny],changevec(:,1),changevec(:,2)))=changevec(:,3);
          
            end;
			% Add to the result, multipled by the value in the stencil
			Matt = Matt + tmpmat*stencil(ix+szx+1,iy+szy+1);
		end;
	end;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%% AUX FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%

function changevec = CorrectBC(ix,iy,nx,ny,BC)
% Function to correct for BC's in retrospect
% changevec has a set of [x y z] rows, where x,y are coordinates in the SM,
% and z is the value to change. The third dimension of changevec goes over
% the "problematic" points and where they "end up". the first one is for
% deleting the original sites that contained the periodic information.

vec=[]; % find border sites affected by BC for this site in the stencil
if(ix)  % add abs(ix)*ny sites   
    vec =  [vec; [repmat((1:abs(ix))'+(ix>0)*(nx-ix),ny,1) reshape(repmat((1:ny),abs(ix),1),ny*abs(ix),1)]];
end;
if(iy)  % add abs(iy)*nx-abs(ix)*abs(iy) sites   
    vec =  [vec; [reshape(repmat((1:nx-abs(ix))+(ix<0)*abs(ix),abs(iy),1),(nx-abs(ix))*abs(iy),1) repmat((1:abs(iy))'+(iy>0)*(ny-iy),nx-abs(ix),1)]];
end;
if(~isempty(vec))
    vecdst = [mod(vec(:,1)+ix-1,nx)+1 mod(vec(:,2)+iy-1,ny)+1];
    %disp([vec vecdst])
    changevec = [sub2ind([nx,ny],vecdst(:,1),vecdst(:,2))  sub2ind([nx,ny],vec(:,1),vec(:,2)) -ones(size(vecdst,1),1)];
	if(BC>0)    % Neumann BC
        [newvec,newparts]=MoveInStencil([vec(:,1)+ix vec(:,2)+iy],[ix iy],1,nx,ny);
        %tmp = [vecdst -ones(size(vecdst,1),1)];
        for ii=1:length(newparts)
            %tmp = cat(3,tmp,[newvec(:,:,ii) repmat(newparts(ii),size(vecdst,1),1)]);
            tmp = [sub2ind([nx,ny],newvec(:,1,ii),newvec(:,2,ii)) sub2ind([nx,ny],vec(:,1),vec(:,2)) repmat(newparts(ii),size(vecdst,1),1)];
            changevec = cat(3,changevec,tmp);
        end;
        %vecdst = [mod(vec(:,1)+ix-1,nx)+1 mod(vec(:,2)+iy-1,ny)+1];
        %changevec = [  sub2ind([nx,ny],vecdst(:,1),vecdst(:,2))  sub2ind([nx,ny],vec(:,1),vec(:,2)) -ones(size(vecdst,1),1)];
    end;
else
	changevec=vec;
end;
   
end

function [newvec,newparts]=MoveInStencil(vec,loc,part,nx,ny)
% Recursive function to slowly bring sites from some site in the stencil 
% to the center 0,0 site, collecting information along the way.
ix=loc(1);
iy=loc(2);

% Which sites still need to be taken in
needchange = ~(prod([mod(vec(:,1)-1,nx)+1 mod(vec(:,2)-1,ny)+1]==vec,2)); 

% Did we reach the end? (all sites in "valid" places = no more corrections)
if(sum(needchange)==0)
    newvec = vec;
    newparts = part;
else
    if(ix&&iy) % Divide into two paths, (on x and y), with part being split
        tmpvec=vec;
        tmpvec(needchange,1)=tmpvec(needchange,1)-sign(ix);
        [newvec1,newparts1]=MoveInStencil(tmpvec,[loc(1)-sign(ix) loc(2)],part*abs(ix)/(abs(ix)+abs(iy)),nx,ny);
        tmpvec=vec;
        tmpvec(needchange,2)=tmpvec(needchange,2)-sign(iy);
        [newvec2,newparts2]=MoveInStencil(tmpvec,[loc(1),loc(2)-sign(iy)],part*abs(iy)/(abs(ix)+abs(iy)),nx,ny);
        newvec = cat(3,newvec1,newvec2);
        newparts = [newparts1;newparts2];
        
    else % Go into either x or y path (which is ever is not 0)
        tmpvec=vec;
        tmpvec(needchange,1)=tmpvec(needchange,1)-sign(ix);
        tmpvec(needchange,2)=tmpvec(needchange,2)-sign(iy);
        [newvec,newparts]=MoveInStencil(tmpvec,[loc(1)-sign(ix),loc(2)-sign(iy)],part,nx,ny);
    end;
end;
    
end