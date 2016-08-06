function [P2Pout,resout] = PrepSol_P2P2(Vs,Ps,Es,InitFuncP2P,varargin)
% prepare a new variable for pde2path(V2) continuation, using given data

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

MaxNewtonIter = 200;
%if(nargin<5)
%    name = 'pp';
%end
if(~iscell(InitFuncP2P))                        % Force into cell structure
    InitFuncP2P={InitFuncP2P};
end;

pp=[];
pp=InitFuncP2P{1}(pp,InitFuncP2P{2:end});       % run init function

previmax = pp.asw.imax;
pp.asw.imax = MaxNewtonIter;

%ChangeRes([st Ps_LL.SpaMat{1}*st/4]
temp1 = [Vs Ps.SpaMat{1}*Vs];
temp2 = repmat([temp1 ; temp1(1,:)],2,1);
pp.u(1:size(temp2(:))) = temp2(:);


[uu,res,iter,gg1,gg2]=nlooppde(pp,pp.u);
pp.u=uu;

if(res>pp.asw.tol)
    warning('Newton iterations (#%d) failed to reach required res of %e, only: %e',iter,pp.asw.tol,res);
else
    disp(iter);
end;

    prevds = pp.sol.ds;
    prevnsteps = pp.asw.nsteps;
    pp.sol.ds = pp.sol.ds/100;
    pp.asw.nsteps = 1;
    
    try
			pp=cont(pp);
	catch
			warning('cont crashed with previous res: %e',res);
	end;
    %pp = cont(pp);
    
    pp.sol.ds = prevds;
    pp.asw.nsteps = prevnsteps;


pp.asw.imax = previmax;

P2Pout = pp;
resout = res;
end
