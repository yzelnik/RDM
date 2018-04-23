function rhs=RightHandSide(Vs,Ps,Es)
% Return right-hand-side of equation (both local and spatial parts)	

Vs=Vs(:,:,1);
if(Es.SpaMatUse)   % Do we use a spatial matrix?
    rhs = (Ps.LocFunc(Vs,Ps,Es) + reshape(Ps.SpaMat*Vs(:),Ps.Nx*Ps.Ny,Ps.VarNum));     
else        % if we don't use the SM (spatial matrix) than use the spatial function directly
    rhs = (Ps.LocFunc(Vs,Ps,Es) + Ps.SpaFunc(Vs,Ps,Es));
end;

end
