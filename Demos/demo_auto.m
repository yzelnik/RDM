%% DEMO: Using AUTO results
%  In this demo we will learn how to import results from AUTO and use them
%  external files: c.POMK, POMK.c, POMK.auto
%  
%  To start, run POMK.auto using AUTO to get data files we will use below

clear all;
clc;

%% Define model and such
Ps = struct('LocFunc',@L_MK,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'a',1,'m',0.5,'Ds',[1 8],'VarNum',2,'Lx',200,'Ly',1,'Nx',1000,'Ny',1);
Es = struct('TsSize',0.1,'TimeDst',200,'SsThresh',1e-6,'NonNeg',1,'LsaThresh',1e-3,'StSmall',0.01,'VarInd',1);
%  We look at the modified Klausmeier model, of 2 variables B,W:
%  dB/dt =  W*B*B - m*B + B_xx
%  dW/dt =  a - W - W*B*B + d * W_xx
%  Parameters are: a, m, d, with values of: (1,0.5,8)

%% Read the bifurication data from auto and plot it

%  Read the bif data for uniform, periodic and snaking (localized) branches
bfunf = ReadAutoBif('b.MK_unf');
bfper = ReadAutoBif('b.MK_per');
bfsnk = ReadAutoBif('b.MK_snk');

clf;
%  Plot all three branches (combine using a cell-array)
plotbf({bfunf,bfper,bfsnk});
xlim([0.96 1.03]);
%  Note, we do not really need to use Ps/Es here

%% Read states from a branch found with AUTO, and look at some of them

st = ReadAutoStates('s.MK_snk',Ps,Es);

%  Plot out 2 of these states (#3 & #15), centering one of them
subplot(1,2,1);
plotst(st,Ps,Es,3)
subplot(1,2,2);
plotst(M_CenterSt(st,Ps,Es),Ps,Es,15)

%% Now, take one such state and simulate dynamics for it

out2 = run2ss(st(:,:,3),Ps,Es,'Es.OlDraw',1,'Ps.a',1.006);
%  Here we pushed a localized state outside of the snaking range
%  Thus, causing a gradual shift to a periodic state

%% Analysis of multiple states for a bifurication diagram
%  In this case, we want to calculate the stability of a bifurication diagram, 
%  using the T_LSA function that performs numerical linear stability analysis

indices = [1 2 8];
updates = {'a','','Lx'};

stab_snk = AnalyzeAutoStates('MK_snk',Ps,Es,@T_LSA,[],indices,updates);
%  The first input gives the basic file name, so the function reads b.MK_snk and s.MK_snk
%  indices & updates give the columns to read in the b.file, and what
%  variables they stand for in the Ps structure (so that 'a' translates to Ps.a)
%  Thus, each state is read from the s.file, while its parameter values are
%  updated from the b.file, when using the function T_LSA

plotbf(stab_snk)


%% Calculating spatial dynamics eigenvalues for the uniform branch

%  We read the bif file for the uniform branch, reading columns of [a,B,W]
bfunf = ReadAutoBif('b.MK_unf',[1 3 4]);

%  Using runpar, we run a test of calculating the eigenvalues per point
%  (skipping most points, using only 1 point every 10 points)
[~,bftry]=runpar([],Ps,Es,'Es.BfPrm',{'a','Vs(1)','Vs(2)'},'Es.BfRange',bfunf(1:10:3000,:),'Es.FuncList',{@T_SpatialDynamicsEigenvalues},'Es.InitFunc',@M_InitUnfSt);
%  Note that using Es.BfPrm & Es.BfRange, we update both Ps.a and Vs at every point

%  Columns 4-7 in bftry contain 4 eigenvalues that we want to look at
%  We ask if all real/imaginary parts of the 4 eigenvalues are zero
sev = (prod(real(bftry(:,4:7)),2)<eps) + (prod(imag(bftry(:,4:7)),2)<eps) + 1;
%  Concat this result with our original bf data
bftot = [bftry sev];

%  Plot shows 3 different regions by eigenvalue data stored in sev.
%  This also corresponds to stability information we know: solid for stable, 
%  dots for unstable to non-uniform perturbations, dashed for generally unstable
plotbf(bftot,Es,'Es.BfPhases',[1 2 3],'Es.BfStyle',['- ';': ';'--']);
axis([0.995 1.035 0.7 1.4])

%% Create full bifurication diagram

indices = [1 2 8];
updates = {'a','','Lx'};

%  Run analysis on AUTO files, similiarly to before
stab_snk = AnalyzeAutoStates('MK_snk',Ps,Es,@T_LSA,[],indices,updates);
%  Here, the states for the periodic solutions are small in size
%  So we want less points per period (200), but we want to check the stability 
%  for a larger system. We look at 10 periods, by choosing Es.MultPeriod=10.
stab_per = AnalyzeAutoStates('MK_per',Ps,Es,@T_LSA,[],indices,updates,'Es.MultPeriod',10,'Ps.Nx',200);
%  Since there is no defined size for the uniform branch in the run we did,
%  we don't use the third variable/column in indices&updates, and change Ps.Lx
stab_unf = AnalyzeAutoStates('MK_unf',Ps,Es,@T_LSA,[],indices(1:2),updates(1:2),'Ps.Lx',50);

%  Plot the three branches together, with the stability info we added
plotbf({stab_unf,stab_per,stab_snk},Es);
axis([0.99 1.03 0.7 1.4])

%% Now try to re-phase the data, and plot the new result
%  Use the stability data with low-resolution (of the 'a' parameter)
%  on the high-resolution data of the original bif files
%  For example, we take stability info from stab_unf and add it as the last column of bfunf
finunf=RephaseBf(bfunf,stab_unf,Es);
finper=RephaseBf(bfper,stab_per,Es);
%  For the snaking branch, since it goes back and forth, the rephasing can
%  get problematic. We therefore use a slower but more accurate method
finsnk=RephaseBf(bfsnk,stab_snk,Es,'Es.RephaseMode',1);

%  Plot the three branches together, with the stability info we added
plotbf({finunf,finper,finsnk},Es);
axis([0.99 1.03 0.7 1.4])

