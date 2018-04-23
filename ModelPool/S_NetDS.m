function MatOut=S_NetDS(~,Ps,Es)
% Spatial Matrix for dispersal (per-site) across a network 
% MatOut=S_NetDS(Vs,Ps,Es)
% Ps.Nx gives the number of sites, Ps.Ds the diffusion rates per agent type
% The netork is either given by Ps.Net (Ps.Nx by Ps.Nx, diagonal set to zeros)
% or a default of all sites connected. 
% Unlike S_NetDL, the sum of each column is normalized to 2,
% so that the more connected sites get more incoming dispersal
% Normalization means that for a chain network the dispersal is the same as
% for a 1D system with diffusion and  Ps.Lx=Ps.Nx (and also for the S_NetDL)

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Pre caclculate spatial matrix, for future use
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
   len=Ps.Nx;
   if(~isfield(Ps,'Net') || length(Ps.Net)<=1) 
       Ps.Net = ones(len)-diag(ones(len,1)); % Default fully-connected
   end;
   Net = logical(Ps.Net-diag(diag(Ps.Net))); % Setting diagonal to zero
   Net = Net./repmat(max(1,sum(Net,1)),len,1); % Normalizing each row to 1
   Net = 2*sparse(Net-diag(ones(len,1))); % main diagonal is at (-2)

   for ii=1:Ps.VarNum	% Put in derivative sub-matrix in a block-diagonal fashion
        SM((ii-1)*len+(1:len),(ii-1)*len+(1:len)) = Net*Ps.Ds(ii);
   end;
end

