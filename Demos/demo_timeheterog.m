%% DEMO: Temporal heterogeneity 
%  This demo will demonstrate how to use temporal heterogeneity in RDM
clear all;
clc;

%% Define model and such
Ps = struct('LocFunc',@L_MK,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'a',1.01,'m',0.5,'Ds',[1 8],'VarNum',2,'Lx',200,'Ly',1,'Nx',400,'Ny',1);
Es = struct('TsSize',0.5,'TimeDst',200,'SsThresh',1e-6,'NonNeg',1,'LsaThresh',1e-3,'StSmall',0.01,'VarInd',1,'StAxis',[0 2],'OlDraw',1);
%  We look at the modified Klausmeier model, of 2 variables B,W:
%  dB/dt =  W*B*B - m*B + B_xx
%  dW/dt =  a - W - W*B*B + d * W_xx
%  Parameters are: a, m, d, with values of: (1,0.5,8)

% start by getting a localized state
loc = run2ss([1;0],Ps,Es,'Es.InitFunc',@M_InitMixSt);

%% Periodic forcing in time

% define the frames, and then sinusodial forcing
frms   = 0:1:800;
sinval = sin(16*2*pi*(frms)/max(frms))*0.03;

% run frames with the sinusodial forcing of parameter a
out = runframes(loc,Ps,Es,'Es.Frames',frms,'Es.DynPrm','a','Es.DynVal',Ps.a+sinval);
% slowly, more gaps grow & attach to the initial one

% display the simulation and save results for next step
plotwf(out,Ps,Es);
per=out(:,:,end);

%% Repeated singular events
% define a recurring function for cutting down biomass:
Es.RecurFunc=@M_CutVar;
% parameters for this function: [cut down -2, of var 1, 0.08 part of space, randomly] 
Es.ModPrm=[-2 1 0.08 -1]; 

% repeat disturbance every 40 frames, 50 times:
recur = (1:50)*40;
out2 = runframes(per,Ps,Es,'Ps.a',1.015,'Es.Frames',frms,'Es.RecurFrames',recur);
% dynamics of repeated disturbance and system response & recovery

%% Integrating with white noise
%  for using white noise and similar types of temporal heteroneity,
%  a different integration scheme needs to be used

out3 = run2ss(per,Ps,Es,'Ps.NoiseType',[0.2;0.5;0],'Ps.IntegFunc',@I_NoiseDFM,'Es.TsSize',0.002);
%  demographic noise induced on the system

