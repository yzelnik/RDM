function product = T_ProdConRes(Vs,Ps,Es)
% Calculate the overall ecossytem productivity for the Consumers-Resources model
% avgvals = T_ProdConRes(Vs,Ps,Es)

% Initilization
    nlen = Ps.VarNum-Ps.ResNum;
    rlen = Ps.ResNum;
    syssz = Ps.Nx*Ps.Ny;
    N=Vs(:,1:nlen); 
    R=Vs(:,nlen+(1:rlen)); 
    
product = sum(sum(Ps.e.*Ps.c.*repmat(N,1,rlen).*reshape(repmat(R',1,nlen)',syssz,nlen*rlen)))/syssz;

end