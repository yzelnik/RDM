function MatOut=S_NLL(Vs,Ps,Es,varargin)
% Spatial Matrix of (negative) Lejeune-Lefever model
% MatOut=S_NLL(Vs,Ps,Es)
% Ps.Nx (Ps.Ny),Ps.Lx (Ps.Ly) give the number of sites and system size for creating the matrix

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
len=Ps.Nx*Ps.Ny;

B = Vs(:,1,1);


if(Es.Jcob==0)
% Model equations
MatOut = - ( 0.5*(Ps.L^2 - B).*(DervSM(2,Ps,Es)*B) - 0.125*B.*(DervSM(4,Ps,Es)*B) );


else
% Jacobian of equations
base2 = DervSM(2,Ps,Es);			% Get 2nd derivative spatial matrix
base4 = DervSM(4,Ps,Es);			% Get 4th derivative spatial matrix

MatOut = 0.5*base2*Ps.L^2;			% Start with the normal diffusion part

part1 = sparse(B*ones(1,len)).*base2;		% B0 * Derv2 (dB)
part2 = sparse(diag(base2*B));			% dB * Derv2 (B0)
MatOut = MatOut - 0.5*(part1 + part2);		% Add the non-linear diffusion part

part1 = sparse(B*ones(1,len)).*base4;		% B0 * Derv4 (dB)
part2 = sparse(diag(base4*B));			% dB * Derv4 (B0)
MatOut = - (MatOut - 0.125*(part1 + part2));	% Add the non-linear 4th derivative part

end

