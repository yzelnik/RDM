function Mat=S_Diffusion(Vs,Ps,Es)
% we return a diffusion-operator matrix
% such that if we call Mat that we return M, and Ps.Ds we call d, we have:
% d*u_xx = M*u

% We begin by building a second-derivative (operator) matrix
% such a second-derivative matrix in 1D looks much like:
% -2 1  0 0 . . 
%  1 -2 1 0 . . 
%  0 1 -2 1 . . 
%  . .  . . . . 

% To create a mostly-zeros matrix, with a few diagonals, we can use:
tmp  = ones(Ps.Nx,1);
base = -2*diag(tmp)+circshift(diag(tmp),1) + circshift(diag(tmp),-1);

% We calculate the spatial resolution
res = Ps.Lx/Ps.Nx;
% Finally, we multiply the matrix by the diffusion coefficient (Ps.Ds),
% and divide by the square of the resolution, to get our diffusion operator
Mat = base * Ps.Ds /res^2;

end