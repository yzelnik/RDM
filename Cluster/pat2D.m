function pat2D(pp)

addpath(genpath('/fastspace/users/kinast/NFW'));

Es=struct...
('Tstep',0.05,'Tdest',200,'Nsteps',1e2,'Jcob',0,'posflag',1,'LSAthresh',0,'STsmall',1e-4,...
'Vind',[1 2],'OLdraw',1,'St1Color',circshift(hsv(4),-1)*0.66,...
'BfColor',[0 0 0 ; 1 0 0; 0 0.75 0; 0 0 1; 0 0.75 0.75; 1 0 1],'BfStyle',['-- ';'-  '],...
'TestFuncs',{{@T_L2Norm,@T_CalcWL,@T_CountRegions,@T_MinMax}},'TestOuts',[2,3,1,4],'Frames',100,...
'InitParms',[5 1],'InitByODE',1);
Es.Frames = [1 5 10 50 100 500 1000 2999];
%Es.Frames = 10;

Ps_DC=struct...
('LocFunc',@L_DC,'SpaFunc',@S_RD,'P',170,'LambdaB',0.032,'EE',1.5,'KB',1,'MB',1.2,'PhiB',1,'beta',1,...
'BZero',0.05,'nn',1.0,'mm',1.0,'LambdaCW',0.035,'LambdaCH',0.01,'KC',0.003,'MC',0.2,'PhiC',20,'A',10,'QC',0.0006,'fC',0.1,'QB',0.05,...
'fB',1.0,'N',4,'R',0.95,'GammaB',30,'GammaCW',0.1,'GammaCH',0.02,'Vnum',4,...
'Ds',[6.25*10^(-4) 6.25*10^(-3) 6.25*10^(-2) 5*10^(0)],'Lx',20,'Ly',20,'Nx',512,'Ny',1);
Ps_DCn=ConvPars_DC([],Ps_DC,[]);
Ps_DCn2=Ps_DCn;
nn=256;Ps_DCn2.Nx=nn;Ps_DCn2.Ny=nn;
ll=500;Ps_DCn2.Lx=ll;Ps_DCn2.Ly=ll;

u2 = InitUnfSt(1,Ps_DCn2,Es,'Es.StNoise',1e-1,'Ps.p',pp);
tic;
w2 = GetSnapshots(u2,Ps_DCn2,Es,@IPS_AKK,'Es.Tstep',0.1,'Ps.p',pp,'Es.Tdest',3000);
toc;

save(['pat2D_' num2str(pp*100)],'w2');

end

