%% Using 1D and 2D systems, getting different spatial patterns

Ps=struct('LocFunc',@L_GS,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'f',0.06,'k',0.06,'Ds',[1 10],'VarNum',2,'Lx',200,'Ly',1,'Nx',200,'Ny',1);
Es=struct('TsSize',0.2,'TimeDst',200,'OdeInit',1,'SsThresh',1e-6,'NonNeg',1,'LsaThresh',1e-3,'StSmall',0.01,'VarInd',1,'StAxis',[0 1]);

%% getting a few peaks in 1D

%% getting a localized pattern in 1D

%% getting stripes in 2D

%% getting hexagonal-gaps in 2D

%% getting a pattern-gradient

%% using a mask for initial conditions