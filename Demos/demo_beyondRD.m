%% DEMO: Going beyond Reaction-Diffusion
%  In this demo we will look at different options of spatial functions beyond Reaction-Diffusion

clear all;
clc;
%%  We start with linear derivatives, but not Reaction-Diffusion

%  For example, we look at the Swift-Hohenberg model, that can be defined as:
%  dU/dt = lambda*u + a2*u^2 + a3*u^3  + d2* D^2(u) + d4* D^4(u)
%  Where D^2 and D^4 are the second and fourth spatial derivatives
%  Parameters are: lambda,a2,a3,d2,d4, with values of: (0,2,-1,-2,-1)

Ps1 = struct('LocFunc',@L_SH,'SpaFunc',@S_LD,'IntegFunc',@I_PSRD,'lambda',-0.9,'a2',2,'a3',-1,'Ds',[0 -2 0 -1],'VarNum',1,'Lx',50,'Ly',50,'Nx',500,'Ny',1);
Es  = struct('TsSize',0.1,'TimeDst',200,'SsThresh',1e-6,'LsaThresh',1e-3,'StSmall',0.01,'VarInd',1,'StAxis',[-1 2.5]);

%  Looking for a non-trivial uniform solution
getode(1,Ps1,Es)

%% Find a non-uniform solution in 1D
out1 = run2ss(1,Ps1,Es,'Es.OlDraw',1);

%% Run the same thing, but with different boundary conditions
out2 = run2ss(1,Ps1,Es,'Es.OlDraw',1,'Ps.IntegFunc',@I_FDSIMP,'Ps.Bc',-1);
%  To do this we cannot use I_PSRD as the integrating function. Note that only
%  now is the choice of S_LD as the spatial function becomes meaningful (i.e. it is now used)

%% Run a 2D simulation for the SH model
out3 = run2ss(1,Ps1,Es,'Es.OlDraw',1,'Ps.Nx',100,'Ps.Ny',100,'Es.St2Interp',1);
%  When both Ps.Nx>1 and Ps.Ny>1, a 2D configuration is implicitly assumed

%% We can also look at a Reaction-Diffusion-Advection model
%  We use the GS model, see the demo_mainfuncs for more explanations

Ps2 = struct('LocFunc',@L_GS,'SpaFunc',@S_RD,'IntegFunc',@I_FDCN,'f',0.06,'k',0.06,'Ds',[1 10],'VarNum',2,'Lx',200,'Ly',1,'Nx',200,'Ny',1);

%  With just normal RD, we get a periodic state with a few peaks
rd  = run2ss(1,Ps2,Es,'Es.OlDraw',1,'Es.StAxis',[0 1]);

%% To add advection, we change the spatial function to S_LD, and so we need 
%  to change the derivative coefficients. The first 2 values of Ps.Ds are for 
%  the first derivative (we choose the have advection for the second variable),
%  while the next two are for the second derivatives (which we keep the same)
rda = run2ss(1,Ps2,Es,'Es.OlDraw',1,'Ps.SpaFunc',@S_LD,'Ps.Ds',[0 0.25 1 10],'Es.StAxis',[0 1]);

%% Non-linear derivatives
%  We can also use non-linear spatial derivatives, which limits our options
%  of integration functions, among other things. (No spectral/implicit methods)

%% The Burger's equation 
%  We use S_BE for a spatial function (with no local terms) to get the Burger's equation
Ps3 = struct('LocFunc',@L_Null,'SpaFunc',@S_BE,'IntegFunc',@I_FDE,'f',0.06,'k',0.06,'Ds',[4 1],'VarNum',1,'Lx',200,'Ly',1,'Nx',1000,'Ny',1);

try1 = run2ss([1;0],Ps3,Es,'Es.TsSize',2*1e-3,'Es.OlDraw',1,'Es.InitFunc',@M_InitMixSt,'Es.InitPrm',0.2);
%  Note that the S_BE function returns a vector the size of Vs, not a spatial matrix

%% Integral terms
%  Using integral terms is typically quite heavy on computing time.


