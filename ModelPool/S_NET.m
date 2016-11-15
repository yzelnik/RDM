function MatOut=S_NET(~,Ps,Es)
% Spatial Matrix for a set of sites connected by a network 
% MatOut=S_NET(Vs,Ps,Es)
% Ps.Nx gives the number of sites, Ps.Ds the diffusion rates per agent type
% The netork is either given by Ps.Net (Ps.Nx by Ps.Nx, diagonal set to zeros)
% or use a default of all-connected. Sum of each column is normalized to 1.

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Pre caclculate spatial matrix, for future use
    MatOut = Ps;
    MatOut.SpaMat = CalcSM(Ps,Es);

else  	 % Normal run
   if(~isfield(Ps,'SpaMat'))    % Caclculate spatial matrix if needed
        Ps.SpaMat = CalcSM(Ps,Es);
   end;
   
   MatOut = Ps.SpaMat;
end

end


%------------ AUX funcs ----------------------------------------------------

function SM=CalcSM(Ps,Es)
   len=Ps.Nx;
   if(~isfield(Ps,'Net') || length(Ps.Net)<=1)
       Ps.Net = ones(len)-diag(ones(len,1));
   end;
   Net = Ps.Net./repmat(sum(Ps.Net,1),len,1); % Normalizing each row to one
   Net = sparse(Net-diag(ones(len,1)));
  
   for ii=1:Ps.VarNum	% Put in derivative sub-matrix in a block-diagonal fashion
        SM((ii-1)*len+(1:len),(ii-1)*len+(1:len)) = Net*Ps.Ds(ii);
   end;

end

