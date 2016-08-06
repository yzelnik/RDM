function MatOut=S_RDWW(Vs,Ps,Es)
% Spatial Matrix of Reaction Diffusion system with wind damage/protection
% MatOut=S_RDWW(Vs,Ps,Es)
% Ps.Nx (Ps.Ny),Ps.Lx (Ps.Ly) give the number of sites and system size for creating the matrix
% The array Ps.Ds is assumed to contain the diffusion rates of the different reactants 
% Ps.xi controls the force of wind protection, Ps.rad is its radius
% B motality is modified to be M/(1+xi*int(B)/(S*K)) where int is an integral done over radius rad
% K(or kappa) and M(or mu) are capacity rate and mortality, assumed to exist or equal to 1, and S is the integration area

if(~isfield(Es,'fmod'))
   Es.fmod=0;
end;

if(Es.fmod<0)	% Pre caclculate spatial matrix, for future use
   MatOut = Ps;
   MatOut.SpaMat = CalcSM(Ps,Es);
   MatOut.IntMat = CalcIntMat(Ps,Es);
else		% Normal run
   if(~isfield(Ps,'SpaMat') || ~isfield(Ps,'IntMat') || isempty(Ps.SpaMat) || isempty(Ps.IntMat)) 
   	Ps.SpaMat = CalcSM(Ps,Es);
	Ps.IntMat = CalcIntMat(Ps,Es);
   end;
   if(isfield(Ps,'K'))
	kappa = Ps.K;
   elseif(isfield(Ps,'kappa'))
	kappa = Ps.kappa;
   else
	kappa = 1;
   end;
   if(isfield(Ps,'M'))
	mu = Ps.M;
   elseif(isfield(Ps,'mu'))
	mu = Ps.mu;
   else
	mu = 1;
   end;

   %size(Ps.SpaMat),size(Vs(:))
   MatOut = Ps.SpaMat*(Vs(:));		% Put in derivatives
   MatOut = reshape(MatOut,size(Vs));

   % Add in the integral, and correct for the old local term by adding in the "-1" term in the end
   MatOut(:,1) = MatOut(:,1) - mu*(1./(1+(Ps.xi./kappa)*Ps.IntMat*Vs(:,1)).^4 -1).*Vs(:,1);

end

end

%------------ AUX funcs ----------------------------------------------------

function SM=CalcSM(Ps,Es)	% This is a normal derivative function
   Es.fmod=-1;
   tempPs = S_RD([],Ps,Es);	% Copy the normal spatial matrix from the RD function
   SM = tempPs.SpaMat;
   Es.fmod=0;
end


function SM=CalcIntMat(Ps,Es)	% Integration matrix
if(Ps.Nx==1)	% This handles only for 1D, should be changed...
   len  = Ps.Ly;
   pnum = Ps.Ny;
elseif(Ps.Ny==1)
   len  = Ps.Lx;
   pnum = Ps.Nx;
end;
rad = round(min(Ps.rad,len/2.5)/len*pnum);

xx  = zeros(pnum,1);
xx(1:rad+1)=1;
xx(end-rad+1:end)=1;

for ii=1:pnum 
   SM(:,ii)=circshift(xx,ii); 
end; 
SM = SM/(rad*2+1);	
end


