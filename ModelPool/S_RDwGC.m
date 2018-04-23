function MatOut=S_RDwGC(Vs,Ps,Es)
% Reaction Diffusion with Global Competition
% MatOut=S_RDwGC(Vs,Ps,Es)
% Ps.Nx (Ps.Ny), Ps.Lx (Ps.Ly) give the number of sites and system size for creating the matrix
% The array Ps.Ds is assumed to contain the diffusion rates of the different reactants 
% Ps.GlobComp is a matrix with the competition coefficients

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Pre calculate spatial matrix, for future use
    MatOut = Ps;
    if(isfield(Es,'SpaMatUse') && Es.SpaMatUse==-1)
        Es.SetupMode=0;
        MatOut.SpaData=S_RD(1,Ps,Es); % do not use proper spatial matrix
    else    
        MatOut.SpaMat = CalcSM(Ps,Es);
    end;
else  	 % Normal run
   if(~isfield(Es,'SpaMatUse')) % Calculate spatial matrix if needed
        Ps.SpaMat = CalcSM(Ps,Es);
   elseif Es.SpaMatUse % Normal use of spatial matrix
       MatOut = Ps.SpaMat;
   else % not returning spatial matrix, often more efficient
       MatOut = reshape(Ps.SpaData*Vs(:),Ps.Nx*Ps.Ny,Ps.VarNum) - repmat(sum(Ps.GlobComp.*repmat(mean(Vs,1),Ps.VarNum,1),2)',Ps.Nx*Ps.Ny,1);
   end;
   
end

end

%------------ AUX funcs ----------------------------------------------------

function SM=CalcSM(Ps,Es)
len=Ps.Nx*Ps.Ny;
SM = sparse(len*Ps.VarNum,len*Ps.VarNum);	% Build the block-diagonal matrix
dsm2 = DervSM(2,Ps,Es);  
   
if(size(Ps.Ds,2)==1) Ps.Ds=Ps.Ds'; end; % Make sure Ds is in a row shape
   
for ii=1:Ps.VarNum				% Put in derivative sub-matrix in a block-diagonal fashion
	SM((ii-1)*len+(1:len),(ii-1)*len+(1:len)) = dsm2.*Ps.Ds(:,ii);
end;
   
% Add global competition 
ons = ones(len)/len; 
for ii=1:Ps.VarNum
	for jj=1:Ps.VarNum
        SM((ii-1)*len+(1:len),(jj-1)*len+(1:len)) = SM((ii-1)*len+(1:len),(jj-1)*len+(1:len)) - Ps.GlobComp(ii,jj).*ons;
	end;
end;

end

