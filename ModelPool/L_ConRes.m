function VsOut=L_ConRes(Vs,Ps,Es)
% Consumers-Resources model - Local terms
% VsOut=L_ConRes(Vs,Ps,Es)
% Variables are: N_1..N_M and R1..R_L, Parameters are: (e,c,m,I,l)
% Where L = Ps.ResNum (default is 1), and M = Ps.VarNum-Ps.ResNum
% Equations: dN_i/dt = N_i*(e_i*c_i*R_j - m_i) + D_i*N_i''
%            dR_j/dt = I - l*R_j - R_j*sum_i(c_i*N_i_) + D_{i+j}*R_j''
% Diffusion is set by Ds, where 1 value given means same diffusion for all,
% and 2 values give it for consumers and resources seperately (otherwise, all should be defined)
% For other parameters, they can each be defined flexibly, 
% (either as scalar,vector or matrix) where for e&c the size of the
% vector/matrix on the N/R axis can be either one of them, or both (N before R)

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Make some adjusments, to allow for easier definition of parameters
    if(~isfield(Ps,'ResNum') || isempty(Ps.ResNum) || Ps.ResNum<1)
        Ps.ResNum=1;
    end;
    Ps=FixupParameters(Ps);

    VsOut = Ps;
else
    
    % Initilization
    nlen = Ps.VarNum-Ps.ResNum;
    rlen = Ps.ResNum;
    syssz = Ps.Nx*Ps.Ny;
    N=Vs(:,1:nlen); 
    R=Vs(:,nlen+(1:rlen)); 
    
    if(Es.JacMode==0)      % Model equations
        % Product of c*N*R, where N&R are extended as necessary
        nonlinprod = Ps.c.*repmat(N,1,rlen).*reshape(repmat(R',1,nlen)',syssz,nlen*rlen);
        
        dN = sum(reshape(Ps.e.*nonlinprod,syssz,nlen,rlen),3) - Ps.m.*N;
        dR = Ps.I - Ps.l.*R - squeeze(sum(reshape(nonlinprod,syssz,nlen,rlen),2));
    
        VsOut = [dN,dR];
    else    % Jacobian of equations
        % NEEDS TO BE WRITTEN
        NdN = 0;
        NdR = 0;
        RdN = 0;
        RdR = 0;
    
    % written in a large sparse matrix format
    VsOut = 0;%spdiags([VdV VdU UdV; UdU VdU UdV],[0 syslen -syslen],syslen*2,syslen*2);
    %VsOut = sparse([diag(VdV) diag(VdU) ; diag(UdV) diag(UdU)]);
    end;
end;

end

%------------ AUX funcs ----------------------------------------------------

function Ps=FixupParameters(Ps)

nlen  = Ps.VarNum-Ps.ResNum;
rlen  = Ps.ResNum;
syssz = Ps.Nx*Ps.Ny;
% setup number of diffusion coefficients
if(length(Ps.Ds(:))<3) 
    if(length(Ps.Ds(:))==2)
        Ps.Ds = [repmat(Ps.Ds(1),1,nlen) repmat(Ps.Ds(2),1,rlen)];
    else
        Ps.Ds = repmat(Ps.Ds(1),1,Ps.VarNum);
    end;
end;
% Extend parameters m,I,l in either system size of number of consumer/resources
if((size(Ps.m,1)==1)&&(size(Ps.m,2)>1))
    Ps.m=repmat(Ps.m,syssz,1);
elseif((size(Ps.m,1)>1)&&(size(Ps.m,2)==1))
    Ps.m=repmat(Ps.m,1,nlen);
end;
if((size(Ps.I,1)==1)&&(size(Ps.I,2)>1))
    Ps.I=repmat(Ps.I,syssz,1);
elseif((size(Ps.I,1)>1)&&(size(Ps.I,2)==1))
    Ps.I=repmat(Ps.I,1,rlen);
end;
if((size(Ps.l,1)==1)&&(size(Ps.l,2)>1))
    Ps.l=repmat(Ps.l,syssz,1);
elseif((size(Ps.l,1)>1)&&(size(Ps.l,2)==1))
    Ps.l=repmat(Ps.l,1,rlen);
end;


% extend parameters e,c in either/and system size, consumer num, resource num
if((size(Ps.e,1)==1)&&(size(Ps.e,2)>1))
    Ps.e=repmat(Ps.e,syssz,1);
elseif((size(Ps.e,1)>1)&&(size(Ps.e,2)==1))
    Ps.e=repmat(Ps.e,1,nlen*rlen);
end;
if((size(Ps.e,2)>1) && (size(Ps.e,2)<nlen*rlen)) % Check for partial agent-num extension
	if(rlen==nlen)
        error('Since number of consumers & resources is the same, Ps.e needs to be fully specified.'); 
	elseif(size(Ps.e,2)==nlen)
        Ps.e=repmat(Ps.e,1,rlen);
    elseif(size(Ps.e,2)==rlen)
        Ps.e=reshape(repmat(Ps.e',1,nlen)',syssz,nlen*rlen);
    else
        error('Ps.e appears to have a bad structure in the N/R number.');
    end;
end

if((size(Ps.c,1)==1)&&(size(Ps.c,2)>1))
    Ps.c=repmat(Ps.c,syssz,1);
elseif((size(Ps.c,1)>1)&&(size(Ps.c,2)==1))
    Ps.c=repmat(Ps.c,1,nlen*rlen);
end;
if((size(Ps.c,2)>1) && (size(Ps.c,2)<nlen*rlen)) % Check for partial agent-num extension
	if(rlen==nlen)
        error('Since number of consumers & resources is the same, Ps.c needs to be fully specified.'); 
	elseif(size(Ps.c,2)==nlen)
        Ps.c=repmat(Ps.c,1,rlen);
    elseif(size(Ps.c,2)==rlen)
        Ps.c=reshape(repmat(Ps.c',1,nlen)',syssz,nlen*rlen);
    else
        error('Ps.c appears to have a bad structure in the N/R number.');
    end;
end;

end