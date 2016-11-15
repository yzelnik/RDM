%% Using RDM on a cluster
%  The main aim to use runpar using the runclust.sh shell script

%  In general, a good idea of how to do this is to do the following:
%  1) run a single run-func (not using runpar) that you want to do
%  2) now use runpar over a small parameter range, in a similar way
%  3) now use runclust.sh but run it on your local computer
%  4) move the relevant files to the cluster, and run it there.

%% (1) run a single scenario
Ps=struct('LocFunc',@L_GS,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'f',0.06,'k',0.06,'Ds',[1 10],'VarNum',2,'Lx',200,'Ly',1,'Nx',200,'Ny',1);
Es=struct('TsSize',0.2,'TimeDst',200,'OdeInit',1,'SsThresh',1e-6,'NonNeg',1,'LsaThresh',1e-3,'StSmall',0.01,'VarInd',1,'StAxis',[0 1]);

out1 = runflow(1,Ps,Es,'Es.OlDraw',1);

%% (2) run this scenario over a small parameter range
out2 = runpar(1,Ps,Es);

%% (3) build a script file, or some other thing, and run runclust.sh

%% (4) run things on the cluster

