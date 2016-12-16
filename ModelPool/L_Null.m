function VsOut=L_Null(~,Ps,Es)
% No local terms
% VsOut=L_Null(Vs,Ps,Es)

len = Ps.Nx*Ps.Ny;
if(~isfield(Es,'JacMode') || (Es.JacMode==0))	
    VsOut = sparse(len,Ps.VarNum);
else
    VsOut = sparse(len*Ps.VarNum,len*Ps.VarNum);
end
