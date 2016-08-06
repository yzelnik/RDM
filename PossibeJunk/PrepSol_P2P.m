function P2Pout = PrepSol_P2P(Vs,Ps,Es,InitFuncP2P,varargin)
% prepare a new variable for pde2path continuation, using given data

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

MaxNewtonIter = 500;
%if(nargin<5)
%    name = 'pp';
%end
if(~iscell(InitFuncP2P))    % Force into cell structure
    InitFuncP2P={InitFuncP2P};
end;

pp=[];
pp=InitFuncP2P{1}(pp,InitFuncP2P{2:end});
previmax = pp.imax;
pp.imax = MaxNewtonIter;


%ChangeRes([st Ps_LL.SpaMat{1}*st/4]
temp1 = [Vs Ps.SpaMat{1}*Vs];
temp2 = repmat([temp1 ; temp1(1,:)],2,1);
pp.u = temp2(:);


[uu,res,iter,gg1,gg2]=nlooppde(pp,pp.u,pp.lam-0.001);
pp.u=uu;

if(res>pp.tol)
    warning('Newton iterations (#%d) failed to reach required res of %e, only: %e',iter,pp.tol,res);
end;

    prevds = pp.ds;
    prevnsteps = pp.nsteps;
    pp.ds = pp.ds/100;
    pp.nsteps = 1;
    
    try
			pp=cont(pp);
	catch
			warning('cont failed with lambda %f and previoud res: %e',pp.lam,res);
	end;
    %pp = cont(pp);
    
    pp.ds = prevds;
    pp.nsteps = prevnsteps;


pp.imax = previmax;

P2Pout = pp;

end
