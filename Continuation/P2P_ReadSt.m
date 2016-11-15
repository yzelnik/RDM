function [Vs,pp]=P2P_ReadSt(state,Ps,Es,varargin)
% Read P2P state (either a file name or a p2p structure)
% Vs=P2P_ReadSt(state,Ps,Es)

% Update online if necessary
[~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});

%ffname=[pre '/' fname '.mat'];
if(isstruct(state))
    pp=state;
else    
    s = load([state '.mat'],'p'); 
    pp = s.p; 
end;


lens = max(pp.mesh.geo([2 4],:)')-min(pp.mesh.geo([2 4],:)');   % system size

if(Ps.Ny==1)                % fix up the "correct=dummy" system size for 1D case
    Ps.Ly = lens(2);
end;
if(Ps.Nx==1)
    Ps.Lx = lens(2);
end;



if(Ps.Lx~=lens(1))          % Turn system (90 deg) if there's a mis-match
    tmp = Ps.Lx;
    Ps.Lx = Ps.Ly;
    Ps.Ly = tmp; 
    turn = 1;
else
    turn = 0;
end;


if(Ps.Lx~=lens(1) || Ps.Ly~=lens(2))    % If there's still no match, quit
    error('System size in Ps (%.2f x %.2f) does not match p2p one (%.2f x %.2f)',Ps.Lx,Ps.Ly,lens(1),lens(2));
end;

for ii=1:pp.nc.neq                     % Interpolate mesh to grid
    uxy(:,:,ii) = tri2grid(pp.mesh.p,pp.mesh.t,pp.u((1:pp.np)+(ii-1)*pp.np),(1:Ps.Nx)/Ps.Nx*Ps.Lx-Ps.Lx/2,(1:Ps.Ny)/Ps.Ny*Ps.Ly-Ps.Ly/2);
end;
%size(uxy)
if(~turn)                                % Turn 90 deg back, if turned before
    uxy = permute(uxy,[2 1 3]);
    %disp(111);
end;

Vs=reshape(uxy,Ps.Nx*Ps.Ny,pp.nc.neq); % Set to the correct form

if isfield(Es,'TransDef')               % Transform to a new variable set if needed
    Ps.VarNum=pp.nc.neq;
    Vs = TransformVar(Vs,Ps,Es,Es.TransDef);
end;


end
