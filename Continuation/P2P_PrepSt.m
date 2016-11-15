function [P2Pout,resout] = P2P_PrepSt(Vs,Ps,Es,InitFuncP2P,varargin)
% prepare a new variable for pde2path(V2) continuation, using given data

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
MaxNewtonIter = 1000;    % Arbitary limit for newton loop

tic;

%%% Variable transformation
if(Ps.VarNum==1) && (~isfield(Es,'TransDef'))
% If no transformation is defined, but we have only 1 variable, 
% then assume a SH-like model, so transform it to a 2-variable system (u,v) where u=v''
    Es.TransDef = [1 3]; 
end;

if isfield(Es,'TransDef')   % Transform to a new variable set if needed
    Vs = TransformVar(Vs,Ps,Es,Es.TransDef);
    Ps.VarNum = length(Es.TransDef);
end;

%figure(55);
%plotst(Vs(:,2),Ps,Es,1,'Ps.VarNum',1);
%pause;

%%% Deal with p2p init function
if(~iscell(InitFuncP2P))                        % Force into cell structure
    InitFuncP2P={InitFuncP2P};
end;
%pp=[];
%pp=InitFuncP2P{1}(pp,InitFuncP2P{2:end});       % run p2p init function
pp=InitFuncP2P{1}(InitFuncP2P{2:end});       % run p2p init function
toc;
disp('done with init');

%%% Prepare the appropiate Ps for the p2p system
Ps_p2p = Ps;  
% Figure out (in retrospect) what the system size and number of points is
% in the mesh (assuming a simple-rectangle configuration)
lens = max(pp.mesh.geo([2 4],:)')-min(pp.mesh.geo([2 4],:)');   % system size


xplusy = (pp.nu-pp.mesh.nt)/pp.nc.neq-1;       % Find number of points in grid
xtimey = pp.mesh.nt/pp.nc.neq;                 % x=a+b ,  y=a*b
bigsz = xplusy/2+sqrt(xplusy^2-4*xtimey)/2;     % sz1 = (x+sqrt(x^2+4y))/2
smallsz = xplusy-bigsz;                         % sz2 = x-sz1

Ps_p2p.Lx = lens(1);
Ps_p2p.Ly = lens(2);
if(Ps_p2p.Lx>Ps_p2p.Ly)   % Assume that the bigger x/y side has more points
    Ps_p2p.Nx = bigsz;
    Ps_p2p.Ny = smallsz;
else
    Ps_p2p.Nx = smallsz;
    Ps_p2p.Ny = bigsz;
end;

if(Ps_p2p.Lx~=Ps.Lx)        % Turn system (90 deg) if there's a mis-match
    tmp = Ps_p2p.Lx;
    Ps_p2p.Lx = Ps_p2p.Ly;
    Ps_p2p.Ly = tmp; 
    tmp = Ps_p2p.Nx;
    Ps_p2p.Nx = Ps_p2p.Ny;
    Ps_p2p.Ny = tmp; 
end;

%if(Ps_p2p.Ny==2)       % fix up the inherent mismatch in p2p of 2 points for width of 1
%    Ps_p2p.Ny = 1;
%end;
%if(Ps_p2p.Nx==2)       % fix up the inherent mismatch in p2p of 2 points for width of 1
%    Ps_p2p.Nx = 1;
%end;

if(Ps.Ny==1)            % fix up the "correct=dummy" system size for 1D
    Ps.Ly = Ps_p2p.Ly;
end;
if(Ps.Nx==1)
    Ps.Lx = Ps_p2p.Lx;
end;

if(Ps_p2p.Lx~=Ps.Lx || Ps_p2p.Ly~=Ps.Ly)    % If there's still no match, quit
    error('System size in Ps (%.2f x %.2f) does not match p2p one (%.2f x %.2f)',Ps.Lx,Ps.Ly,Ps_p2p.Lx,Ps_p2p.Ly);
end;

% If number of points don't match, interpolate to the new system
if(Ps_p2p.Nx~=Ps.Nx || Ps_p2p.Ny~=Ps.Ny)    
    Vs = ChangeRes(Vs,Ps,Es,Ps_p2p);
    size(Vs)
    disp([Ps_p2p.Nx Ps.Nx  Ps_p2p.Ny  Ps.Ny]);
    disp('interp!');
end;

% Go over each variable, wrap it with the correct edges, and put it in 
for ii=1:Ps_p2p.VarNum
    %size(reshape(Vs(:,ii),Ps_p2p.Nx,Ps_p2p.Ny))
    %size(Vs(1:Ps_p2p.Ny:end,ii))
    %size(Vs(1:Ps_p2p.Ny,ii)')
    %wrap = [reshape(Vs(:,ii),Ps_p2p.Nx,Ps_p2p.Ny) Vs(1:Ps_p2p.Ny:end,ii); Vs(1:Ps_p2p.Ny,ii)' Vs(1,ii)];
    wrap = [reshape(Vs(:,ii),Ps_p2p.Nx,Ps_p2p.Ny) Vs(1:Ps_p2p.Nx,ii); Vs(1:Ps_p2p.Nx:end,ii)'  Vs(1,ii)];
%figure(555);
%    imagesc(wrap);
%    pause;
    p2pu(:,ii)=reshape(flipdim(wrap,2),length(wrap(:)),1); 
end;



pp.u(1:length(p2pu(:))) = p2pu(:);  % Load into the p2p variable

% tri2grid - use??
% % rap4 = reshape(flipdim([reshape(ooo3,20,20,2) ooo3(1:20,:); ooo3(1:21:end,:)' ooo3(1,:)],2),21*21,1,2);
% tmp3 = reshape(vs3,100,100,2);
% ooo3=reshape(tmp3(1:5:end,1:5:end,:),20*20,2);
%temp1 = [Vs];% Ps.SpaMat{1}*Vs];
%temp2 = repmat([temp1 ; temp1(1,:)],2,1);
%plotsol(pp,13,1,2);
%pause;
%plotsol(pp,13,2,2);
%pause;
toc;
tic;
previmax = pp.nc.imax;     % prepare newton loop and run it
pp.nc.imax = MaxNewtonIter;

[uu,res,iter,gg1,gg2]=nloop(pp,pp.u);
pp.u=uu;
toc;
%res=0; iter=1;
if(res>pp.nc.tol)
    warning('Newton iterations (#%d) failed to reach required res of %e, only: %e',iter,pp.nc.tol,res);
else
    disp(iter)
end;

    prevds = pp.sol.ds;
    prevnsteps = pp.nc.nsteps;
    pp.sol.ds = pp.sol.ds/100;
    pp.nc.nsteps = 1;
    
    try         % Just in case cont crashes
			pp=cont(pp);       
	catch
			warning('cont crashed with previous res: %e',res);
	end;
    %pp = cont(pp);
    
    % return the original values
    pp.sol.ds = prevds;
    pp.nc.nsteps = prevnsteps;
    pp.nc.imax = previmax;


P2Pout = pp;
resout = res;

end
