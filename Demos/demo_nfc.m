%% DEMO: some basic simulations of the NFC model
clear all
%% Define the basic structures for our model (model parameters, external parameters)

Ps=struct('LocFunc',@L_VSG,'SpaFunc',@S_RD,'IntegFunc',@I_FDSIMP, 'P',100,'eta',7,'kappa',0.4,'mu',10.5,'nu',15,'lambda',0.9,'gamma',12,'rho',0.7,'Ds',[1.2 150],'VarNum',2,'Lx',150,'Ly',1,'Nx',600,'Ny',1,'Bc',1);
Es=struct('TsSize',0.005,'TimeDst',100,'SsThresh',1e-6,'NonNeg',1,'StSmall',0.01,'VarInd',1,'StAxis',[0 0.3]);

%% simple run, create 8 gaps in 1D

gaps = run2ss([1;0],Ps,Es,'Es.OlDraw',1,'Es.InitFunc',@M_InitPerSt,'Es.InitPrm',[8 0]);

%% simple run, create an invading vegetation front

frnt = run2ss([1;0],Ps,Es,'Es.OlDraw',1,'Es.InitFunc',@M_InitMixSt,'Ps.P',105);

%% get a 2D hexagonal-gap pattern

hex = run2ss([1;0],Ps,Es,'Es.OlDraw',1,'Es.InitFunc',@M_InitPerSt,'Es.InitPrm',[4 -1],'Ps.Lx',40,'Ps.Ly',30,'Ps.Nx',80,'Ps.Ny',60);

%% change conditions (P) periodically, leading to a gradual regime shift

% define the frames, and then sinusodial forcing
frms   = (0:1:240)/2; amp=12; per=8;
rainval = 100+cos(per*2*pi*(frms)/max(frms))*amp;

rng(22); % randomization seed, will affect initial spot locations
% run frames with the sinusodial forcing of percipitation parameter P
spots = runframes([0;1],Ps,Es,'Es.Frames',frms,'Es.DynPrm','P','Es.DynVal',rainval,'Es.OlDraw',2,'Es.InitFunc',@M_InitSpots,'Es.InitPrm',[5 3],'Ps.Lx',50,'Ps.Ly',40,'Ps.Nx',100,'Ps.Ny',80,'Es.TestFunc',@T_CountRegions,'Es.OdeInit',1);
% slowly, more gaps grow & attach to the initial one
% right plot follows the number of gaps over time



