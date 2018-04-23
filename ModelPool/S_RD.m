function MatOut=S_RD(~,Ps,Es)
% Spatial Matrix of Reaction Diffusion system
% MatOut=S_RD(Vs,Ps,Es)
% Ps.Nx (Ps.Ny),Ps.Lx (Ps.Ly) give the number of sites and system size for creating the matrix
% The array Ps.Ds is assumed to contain the diffusion rates of the different reactants 

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Pre calculate spatial matrix, for future use
    MatOut = Ps;
    MatOut.SpaMat = CalcSM(Ps,Es);
else  	 % Normal run
   if(~isfield(Ps,'SpaMat'))    % Calculate spatial matrix if needed
        Ps.SpaMat = CalcSM(Ps,Es);
   end;
   
   MatOut = Ps.SpaMat;
end

end

%------------ AUX funcs ----------------------------------------------------

function SM=CalcSM(Ps,Es)
   len=Ps.Nx*Ps.Ny;
   SM = sparse(len*Ps.VarNum,len*Ps.VarNum);	% Build the block-diagonal matrix
   dsm2 = DervSM(2,Ps,Es);  % second derivative
   
   % Make sure Ds is in a row shape
   if(size(Ps.Ds,2)==1) Ps.Ds=Ps.Ds'; end; 
   
   for ii=1:Ps.VarNum				% Put in derivative sub-matrix in a block-diagonal fashion
        SM((ii-1)*len+(1:len),(ii-1)*len+(1:len)) = dsm2.*Ps.Ds(:,ii);
   end;

end

