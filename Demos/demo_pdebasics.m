%% DEMO: The basics of the pde integration scheme
%  In this demo we will learn the basics of how the PDE integration works
%  external files: L_basics.m, S_Diffusion.m, I_SimpleEuler.m

%  To use the demo, choose each section (by order) with your mouse, 
%  click Cntr-Enter, and see the result, before going to the next section.
clear all;
clc;
%% Define the basic structures for our model (model parameters, external parameters)
%  For Ps, we define 2 functions for the local and spatial parts (L_FKPP,S_RD), 
%  2 model parameters (r,Ds), and 2 parameters for the system/grid size (Lx,Nx)
Ps = struct('LocFunc',@L_FKPP,'SpaFunc',@S_RD,'r',0.5,'Ds',1,'Lx',200,'Nx',400);
%  We also define 2 external-parameters (time-step size, scale of small-values in model)
Es = struct('TsSize',0.01,'StSmall',0.01);

%% Find a steady-state for the non-spatial (ODE) model, starting with a guess u0=0.5
sol = getode(0.5,Ps,Es);
disp(sol)

%% Run a short 1D solution until a steady-state is reached
out1 = run2ss(0.2,Ps,Es);

figure(1); clf;
plotst(out1,Ps,Es);

%% Now, lets create a front, and see its dynamics over time
frnt1 = run2ss([1;0],Ps,Es,'Es.OlDraw',1,'Es.InitFunc',@M_InitMixSt,'Es.StAxis',[-0.1 1.2]);
%  Note: Here we add to Es 3 more parameters: Es.OlDraw for on-line drawing, 
%  Es.InitFunc for initlizing a front, Es.StAxis for keeping a constant y-axis

%% Lets try again, but with a very simple&short local function, @L_basics
Ps2 = Ps; 
Ps2.LocFunc = @L_basics;
%  Where L_basics is a very-short function in this folder. Take a look.

%  Now we run the same front simulation again, and it works just the same
frnt2 = run2ss([1;0],Ps2,Es,'Es.OlDraw',1,'Es.InitFunc',@M_InitMixSt,'Es.StAxis',[-0.1 1.2]);
%  Note that we are using Ps2, which uses the new function

%% But, what about the spatial function and the integration itself?
%  Note that there is a default choice for integration function, that we will overwrite
%  Lets write a very simple Euler-integration scheme
%  Reminder, we assume that: u(t+dt) = u(t) + dt * f(t)
%  Where f(t) is the rhs of our model, or, exactly what's in the L_basics.m function

%% The simple integration function is: @I_SimpleEuler
Ps2.IntegFunc = @I_SimpleEuler; 
frnt2 = run2ss([1;0],Ps2,Es,'Es.OlDraw',1,'Es.InitFunc',@M_InitMixSt,'Es.StAxis',[-0.1 1.2]);

%% And now, all that is left is to define the spatial function/matrix
%  We do this in the function @S_Diffusion. 
%  Look inside the code of @S_Diffusion for more info on that...
Ps2.SpaFunc = @S_Diffusion; 
% We should also define the variable number (due to a technicality)
Ps2.VarNum = 1;
frnt2 = run2ss([1;0],Ps2,Es,'Es.OlDraw',1,'Es.InitFunc',@M_InitMixSt,'Es.StAxis',[-0.1 1.2]);
%  As you might notice, this is much slower, since our code is not optimized.


