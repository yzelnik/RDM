%% DEMO: some basic stuff of the NFC model

%% Define the basic structures for our model (model parameters, external parameters)

Ps=struct('LocFunc',@L_VSG,'SpaFunc',@S_RD,'IntegFunc',@I_FDSIMP, 'P',100,'eta',7,'kappa',0.4,'mu',10.5,'nu',15,'lambda',0.9,'gamma',12,'rho',0.7,'Ds',[1.2 150],'VarNum',2,'Lx',150,'Ly',1,'Nx',600,'Ny',1,'Bc',1);
Es=struct('TsSize',0.005,'TimeDst',100,'SsThresh',1e-6,'NonNeg',1,'StSmall',0.01,'VarInd',1,'StAxis',[0 0.3]);

%% simple run, create 8 gaps in 1D

gaps = run2ss([1;0],Ps,Es,'Es.OlDraw',1,'Es.InitFunc',@M_InitPerSt,'Es.InitPrm',[8 0]);

%% simple run, create an invading vegetation front

frnt = run2ss([1;0],Ps,Es,'Es.OlDraw',1,'Es.InitFunc',@M_InitMixSt,'Ps.P',105);

%% get a 2D hexagonal-gap pattern

hex = run2ss([1;0],Ps,Es,'Es.OlDraw',1,'Es.InitFunc',@M_InitPerSt,'Es.InitPrm',[4 -1],'Ps.Lx',40,'Ps.Ly',30,'Ps.Nx',80,'Ps.Ny',60);

%% a few gaps, using a (1+EB)^1 version now

gaps = run2ss([1;0],Ps,Es,'Ps.LocFunc',@L_VSGE1,'Ps.eta',16,'Es.OlDraw',1,'Es.InitFunc',@M_InitPerSt,'Es.InitPrm',[8 0]);

%% a front, using a (1+EB)^1 version now

frnt = run2ss([1;0],Ps,Es,'Ps.LocFunc',@L_VSGE1,'Ps.P',140,'Es.OlDraw',1,'Es.InitFunc',@M_InitMixSt);
