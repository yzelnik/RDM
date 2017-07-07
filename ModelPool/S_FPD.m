function Out=S_FPD(Vs,Ps,Es)
% Fokker-Planck diffusion
% Out=S_FPD(Vs,Ps,Es)
if(~isfield(Ps,'Nld'))
    Ps.Nld=0;
end;

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Pre caclculate spatial matrix, for future use
    Out = Ps;
    if(~Ps.Nld)
        Out.SpaMat = CalcSM(Ps,Es);
    end;
else  	 % Normal run
    
    if(Ps.Nld)
        len  = Ps.Nx*Ps.Ny;
        SM   = sparse(len*Ps.VarNum,len*Ps.VarNum);	% Build the block-diagonal matrix
        dsm2 = DervSM(2,Ps,Es);  
    
        for ii=1:Ps.VarNum				% Put in derivative sub-matrix in a block-diagonal fashion    
            SM((ii-1)*len+(1:len),(ii-1)*len+(1:len)) = dsm2;
        end; 
        Out = reshape(SM*(Vs(:).*Ps.Ds(:)),len,Ps.VarNum);
    else
        if(~isfield(Ps,'SpaMat'))    % Calculate spatial matrix if needed
            Ps.SpaMat = CalcSM(Ps,Es);
        end;
        Out = Ps.SpaMat;
    end;
end;

end


%------------ AUX funcs ----------------------------------------------------

function SM=CalcSM(Ps,Es)

len  = Ps.Nx*Ps.Ny;
if(size(Ps.Ds,2)==len)
    Ps.Ds=Ps.Ds';
end;
SM   = sparse(len*Ps.VarNum,len*Ps.VarNum);	% Build the block-diagonal matrix
grad = GradSM(1,Ps,Es);
dsm2 = DervSM(2,Ps,Es);  

for ii=1:Ps.VarNum				% Put in derivative sub-matrix in a block-diagonal fashion
    part1 = sparse(Ps.Ds(:,ii)*ones(1,len)).*dsm2;	% D * Derv2 (dU)
    part2 = sparse(diag(dsm2*Ps.Ds(:,ii)));			% dU * Derv2 (D)
    part3 = sparse(len,len);
    for jj=1:length(grad)  % 2* Derv1(D) * Derv1(dU)
       part3=part3 + 2*sparse((grad{jj}*Ps.Ds(:,ii))*ones(1,len)) .*grad{jj};
    end;
    
    SM((ii-1)*len+(1:len),(ii-1)*len+(1:len)) = part1+part2+part3;
end;

  
end

