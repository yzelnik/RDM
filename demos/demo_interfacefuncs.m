%% DEMO: interface functions
%  In this demo we will present the interface functions
%  The names of all these functions start with one of: run, plot, get.
%  The run-funcs are used to perform the main work: run simulations in different modes.
%  The plot-funcs are there to help display the results of runs.
%  The get-funcs are typically used as accessories, to help clarify specific issues.

%  We will present a basic usage of each run-func, by order of increasing complexity
%  The plot/get-funcs will be presented in-between, whereever relevant
clear all;
clc;
%% We start by preparing our 3 basic matlab-variables, using the 1-D Gray-Scott model 
rng(1) % control randomization
Ps=struct('LocFunc',@L_GS,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'f',0.06,'k',0.06,'Ds',[1 10],'VarNum',2,'Lx',200,'Ly',1,'Nx',200,'Ny',1);
Es=struct('TsSize',0.2,'TimeDst',100,'SsThresh',1e-6,'NonNeg',1,'StSmall',0.01,'VarInd',1,'StAxis',[0 1]);
Vs=rand(Ps.Nx,Ps.VarNum);

%% runsim - running a simple simulation pf time-integration

%  This is the simplest simulation type, used to integrate in time until Es.TimeDst
out1 = runsim(Vs,Ps,Es);

%  We can look at the result by using plotst, plotting a state of the system
plotst(out1,Ps,Es);

%% run2ss - running a simulation until steady-state is reached 

%  Here, the runtime is now known before hand, but is determined by the threhshold Es.SsThresh 
out2 = run2ss(out1,Ps,Es,'Es.OlDraw',1);
%  The flag Es.OlDraw indicates we want to see the integration process while it is going on.

%  We use gettest to calculate the wavelength of the solution that we got (size of a single period)
disp(gettest(out2,Ps,Es,@T_CalcWL))

%% runframes - running a simulation while saving some snapshots in the process 

%  Note, we use Es.Frames to give the specific times we want saved
frms = runframes(Vs,Ps,Es,'Es.Frames',1:200);

%  We can look at the result by using plotwf, giving us a waterfall/space-time plot
plotwf(frms,Ps,Es);

%% runnewt - run a Newton-Raphson loop to find a steady-state

%  This is often a fast way to find a steady-state, if we start close enough 
stst = runnewt(out1,Ps,Es);
%  note we use here the output from running runsim

%  plot the input and output of the newton-loop
subplot(1,2,1);
plotst(out1,Ps,Es);
title('before Newton-loop');
subplot(1,2,2);
plotst(stst,Ps,Es);
title('after Newton-loop');


%% runflow - run a few functions consecutively

%  Run a few functions in order: first count the peak-number, then delete a peak, 
%  then integrate in time, then count a peak again. Returns final state and peak numbers.
[out3,stinfo] = runflow(out2,Ps,Es,'Es.FuncList',{@T_CountRegions,@M_DeleteRegions,@runsim,@T_CountRegions},'Es.MergeBfData',1);
%  Note we used the Es.MergeBfData flag to get all the test-results combined

%  Display the test-results of runflow (two different region-counts, before and
%  after removing a region), and also plot out the end result, for context
disp(stinfo);
clf;
plotst(out3,Ps,Es)

%% runfind - iteratively find/optimize some condition

%  find the value of parameter Ps.f for which the average of the first variable is 0.2
[bestst,bestprm] = runfind(Vs,Ps,Es,'Es.FindVal',0.15,'Es.BfPrm','f','Es.FindFunc',@runflow,'Es.FuncList',{@runsim,@T_AvgVal});
%  note that the function we're optimizing is runflow, where we arbitarily running a simulation 
%  for a predefined amount of time, and calculating the average value of variable 1.

%  calculate the average at the original parameters, just for future reference
orgval = gettest(out1,Ps,Es,@T_AvgVal);

%  compare the input and output of runfind
subplot(1,2,1);
plotst(out1,Ps,Es);
title(sprintf('State at origin f=%.3f, avg=%.3f',Ps.f,orgval(1)));
subplot(1,2,2);
plotst(bestst,Ps,Es);
title(sprintf('Found point f=%.3f, with avg=%.3f',bestprm(3),bestprm(1)));

%% runcont - continue a solution, changing a parameter in the process

%  starting with an initial solution, take small steps in parameter f,
%  at each step use a newton-loop, and see how the L2Norm of the solution chages
[finst,branch] = runcont(out1,Ps,Es,'Es.BfPrm','f','Es.BfRange',0.05:0.0005:0.2,'Es.TestFunc',@T_L2Norm);

%  plot a bifurication diagram of a branch that we continued
plotbf(branch,Es);
xlabel('parameter f');
ylabel('L2Norm of the first variable');

%% runpar - run many similar scenarios in parallel (with different parameters)

%  run in parallel a constant-time simulation and count the number of peaks
[~,grid] = runpar(Vs,Ps,Es,'Es.BfPrm',{'k','f'},'Es.BfRange',{[0.02 0.08 7],[0.04 0.18 15]},'Es.FuncList',{@runsim,@T_CountRegions});
%  we run this on a grid of both k&f parameters, where Es.BfRange details the grid used

%%  plot the parameter-space of f/k we get, with color showing how many peaks we have
plotps(grid,Es)

xlabel('k'); ylabel('f'); title('number of peaks');
colorbar;

%% some remaining functions, mostly useful for debugging and similar things:
%  (getode, getrhs, getjac)

%  run getode twice, with different initial guesses for uniform-solution
unf1=getode([0 2],Ps,Es)
unf2=getode([1 1],Ps,Es)

%  run getrhs (and plot it), to get the right-hand-side of the trivial solution
subplot(1,2,1);
plot(getrhs(unf1,Ps,Es))

%  plot our the jacobian of this same state
subplot(1,2,2);
spy(getjac(unf1,Ps,Es))

