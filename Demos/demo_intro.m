%% DEMO: Introduction to RDM
%  This demo will show-case some simple things we can do using RDM
%  It is not meant to explain how this is done "behind the scenes"

%  To use the demo, choose each section (by order) with your mouse, 
%  click Cntr-Enter, and see the result, before going to the next section.
clear all;
clc;
%% We start by defining to basic structures we use in RDM, Ps and Es
Ps=struct('LocFunc',@L_GS,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'f',0.06,'k',0.06,'Ds',[1 10],'VarNum',2,'Lx',200,'Ly',1,'Nx',200,'Ny',1);
Es=struct('TsSize',0.2,'TimeDst',200,'OdeInit',0,'SsThresh',1e-5,'NonNeg',1,'LsaThresh',1e-3,'StSmall',0.01,'VarInd',1,'StAxis',[0 1]);
% Here, Ps defines the parameters of a model (Gray-Scott, in this case)
% While Es defines the external parameters (mostly of a numerical nature)

%% run a short 1D simulation until a steady-state is reached
figure;
out1 = run2ss(1,Ps,Es,'Es.OlDraw',1);
%  On the left panel we see the pattern changing in real time,
%  while on the right panel we follow the change between consecutive frames
%  You can stop the simulation in the middle by clicking on the stop-button

%% Now, run a 2D simulation
out2 = run2ss(1,Ps,Es,'Es.OlDraw',1,'Ps.Nx',80,'Ps.Lx',120,'Ps.Ny',80,'Ps.Ly',120);
%  Note that it takes much longer to get a clear pattern
%  You should probably stop the simulation by pressing the stop-button

%% Lets try to create a front
frnt1 = run2ss([1;0],Ps,Es,'Ps.k',0.033,'Es.OlDraw',1,'Es.InitFunc',@M_InitMixSt,'Ps.Nx',800,'Ps.Lx',800,'Ps.Bc',1,'Ps.IntegFunc',@I_FDSIMP);

%% And in 2D, start with a circle
frnt2 = run2ss([1;0],Ps,Es,'Ps.f',0.03,'Ps.k',0.045,'Es.TsSize',0.5,'Es.OlDraw',1,'Es.InitFunc',@M_InitMixSt,'Ps.Nx',150,'Ps.Lx',300,'Ps.Ny',120,'Ps.Ly',240,'Es.InitPrm',[0.5 0.5 0.1],'Es.StAxis',[0 0.5]);

