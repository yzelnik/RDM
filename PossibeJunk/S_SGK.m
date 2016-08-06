function Out=S_SGK(Vs,Ps,Es,varargin)
% Spatial Matrix of Simplfied Gilad (non-dimensional) model
% VsOut=S_SG(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the spatial terms of the model
% Variables are: B(1),W(2),H(3). 
% Parameters are: P,q,ni,alpha,eta,gamma,rho,f,DW,DH. (1,0.05 3.333 33.333 3.5 16.667  0.95 0.1,100,10000)
% Ps.Nx (Ps.Ny),Ps.Lx (Ps.Ly) give the number of sites and system size for creating the matrix
% The array Ps.Ds is assumed to contain the diffusion rates of the different reactants 

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

len=Ps.Nx*Ps.Ny;

if(length(Ps.SpaMat)<len)
    warning('No Spatial Matrix was defined');
    Ps.SpaMat=DervSM(2,Ps,Es);
end

Out1=Ps.SpaMat*Vs(:,1).*Ps.Ds(1);
Out2=Ps.SpaMat*Vs(:,2).*Ps.Ds(2);
Out3=Ps.SpaMat*(Vs(:,3).^2).*Ps.Ds(3);
Out=[Out1 Out2 Out3];

%for ii=1:Ps.Vnum
%	if(ii<3)
%		MatOut((ii-1)*len+(1:len),(ii-1)*len+(1:len)) = DervSM(2,Ps,Es);
%	else
%		MatOut((ii-1)*len+(1:len),(ii-1)*len+(1:len)) = 2*(DervSM(1,Ps,Es)*Vs(:)*ones(1,len))*DervSM(1,Ps,Es) + 2*(Vs(:)*ones(1,len))*DervSM(2,Ps,Es);
%	end;
%end;
	
end

