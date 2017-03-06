%% DEMO: Spatial heterogeneity 
%  This demo will demonstrate how to use spatial heterogeneity in RDM
clear all;
clc;

Es=struct('TsSize',0.2,'TimeDst',200,'SsThresh',1e-6,'NonNeg',1,'LsaThresh',1e-3,'StSmall',0.01,'VarInd',1);

%% Using a gradient on parameter f in GS model
Ps1=struct('LocFunc',@L_GS,'SpaFunc',@S_RD,'IntegFunc',@I_FDSIMP,'f',0.08,'k',0.06,'Ds',[1 10],'VarNum',2,'Lx',800,'Ly',1,'Nx',800,'Ny',1,'Bc',1);

out = run2ss(0.5,Ps1,Es,'Es.OlDraw',1,'Ps.f',linspace(0.0,0.15,Ps1.Nx)','Es.StAxis',[0 1]);

%% Mixing two gradients in GS model
Ps2=struct('LocFunc',@L_GS,'SpaFunc',@S_RD,'IntegFunc',@I_FDSIMP,'f',0.1,'k',0.06,'Ds',[1 10],'VarNum',2,'Lx',250,'Ly',250,'Nx',100,'Ny',100,'Bc',1);

[gridf,gridk]=meshgrid(linspace(0.0,0.2,Ps2.Nx),linspace(0.0,0.15,Ps2.Nx));
out = run2ss(0.4,Ps2,Es,'Es.OlDraw',1,'Ps.f',gridf(:),'Ps.k',gridk(:));


%% Random heterogeneity (quenched disorder) 
%  changing the carrying capacity (K) in the Logistic model
Ps3=struct('LocFunc',@L_Log,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'r',0.5,'K',1,'Ds',10,'VarNum',1,'Lx',500,'Ly',1,'Nx',200,'Ny',1);

% Run 3 simulations with different levels of heterogeneity
out1 = run2ss(0.5,Ps3,Es,'Ps.K',rand(Ps3.Nx,1)+0.5);
out2 = run2ss(0.5,Ps3,Es,'Ps.K',rand(Ps3.Nx,1)/2+0.75);
out3 = run2ss(0.5,Ps3,Es,'Ps.K',rand(Ps3.Nx,1)/10+0.95);

clf; outs={out1,out2,out3}; % compare the 3 results
for ii=1:3 subplot(1,3,ii); plotst(outs{ii},Ps3,Es,'Es.StAxis',[0 2]); end;

%% Heterogenous diffusion
%  Things get more complicated for heterogeneity in the diffusion itself
%  Here we use a reaction-diffusion setup for Fokker-Planck type of diffusion

Ps4=struct('LocFunc',@L_Log,'SpaFunc',@S_FPD,'IntegFunc',@I_FDSIMP,'r',0.5,'K',1,'Ds',1,'VarNum',1,'Lx',50,'Ly',25,'Nx',100,'Ny',50);
%  The S_FPD assumes Ds has the size of Vs, with the appropiate values of diffusion coefficients 
%  note that the integration function cannot be psedu-spectral, or some other specific type

%  Create a rectangle of higher diffusion values
diffval = ones(Ps4.Nx,Ps4.Ny); 
diffval(11:70,11:30)=2;
diffval = reshape(diffval,Ps4.Nx*Ps4.Ny,1);

out1 = run2ss(1,Ps4,Es,'Es.OlDraw',1,'Ps.Ds',diffval,'Es.TsSize',0.05);

%  plot out the result, showing how the diffusion affects the biomass density
subplot(1,2,1);
plotst(out1,Ps4,Es);
title('biomass density');
subplot(1,2,2);
plotst(diffval,Ps4,Es);
title('diffusion values');
