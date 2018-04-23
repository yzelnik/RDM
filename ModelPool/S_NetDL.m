function MatOut=S_NetDL(~,Ps,Es)
% Spatial Matrix for dispersal (per-link) across a network 
% MatOut=S_NetDL(Vs,Ps,Es)
% Ps.Nx gives the number of sites, Ps.Ds the diffusion rates per agent type
% The netork is either given by Ps.Net (Ps.Nx by Ps.Nx, diagonal set to zeros)
% or a default of all sites connected. 
% Unlike S_NetDR, the sum of each column is not normalized.
% Instead, each link moves agents in rate Ps.Ds so that for a chain network the dispersal 
% is the same as for a 1D system with diffusion and Ps.Lx=Ps.Nx (and also for the S_NetDR)

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
   sumnet = sum(Net,2);
   Net = sparse(Net-diag(sumnet)); % Each column adds up to zero
   
   for ii=1:Ps.VarNum	% Put in derivative sub-matrix in a block-diagonal fashion
        SM((ii-1)*len+(1:len),(ii-1)*len+(1:len)) = Net*Ps.Ds(ii);
   end;
end

