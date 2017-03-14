%% Learning to use more complex runs
%  The main aim to is to learn how to use more complicated scenarios,
%  often when one run-func calls another, that calls another, and such.

%  Define model and such
Ps = struct('LocFunc',@L_MK,'SpaFunc',@S_RD,'IntegFunc',@I_FDSIMP,'a',1.01,'m',0.5,'Ds',[1 8],'VarNum',2,'Lx',200,'Ly',1,'Nx',400,'Ny',1,'Bc',0);
Es = struct('TsSize',0.5,'TimeDst',200,'SsThresh',1e-6,'NonNeg',1,'LsaThresh',1e-3,'StSmall',0.01,'VarInd',1,'StAxis',[0 2],'OlDraw',1);
%  We look at the modified Klausmeier model, of 2 variables B,W:
%  dB/dt =  W*B*B - m*B + B_xx
%  dW/dt =  a - W - W*B*B + d * W_xx
%  Parameters are: a, m, d, with values of: (1,0.5,8)

%% Using runcont which calls runfind

[aba,bla]=runfind([1;0],Ps,Es,'Es.BfPrm','a','Es.TimeMax',200,'Es.OlDraw',0,'Es.FindVal',[1 1e-5 NaN],'Ps.a',1.04,'Ps.Ds',[1 12],'Es.InitFunc',@M_InitMixSt,'Es.InitByOde',1,'Es.SsThresh',1e-5,'Es.TestFunc',@T_CountRegions,'Es.NoWarning',1);
%%
out=run2ss([1;0],Ps,Es,'Ps.a',1.01,'Es.OlDraw',1,'Ps.Ds',[1 8.25],'Es.InitFunc',@M_InitMixSt,'Es.InitByOde',1,'Es.SsThresh',1e-5,'Es.TestFunc',@T_CountRegions);

%% Using runpar which calls runflow which calls run2ss
tic;
[aa,bb]=runcont([1;0],Ps,Es,'Es.BfRange',{[1.04 1.00 51]},'Es.TimeMax',200,'Es.FindOnlyBf',1,'Es.BfOut','MKFrontEndsHere.txt','Es.BfPrm',{'a','Ps.Ds(2)'},'Es.OlDraw',0,'Es.FindVal',[1 1e-5 NaN],'Ps.Ds',[1 13],'Es.InitFunc',@M_InitMixSt,'Es.InitByOde',1,'Es.SsThresh',1e-5,'Es.TestFunc',@T_CountRegions,'Es.SsFunc',@runfind,'Es.NoWarning',1);
toc;
%%

