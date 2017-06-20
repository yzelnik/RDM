function [KMinmax,LambdaMinmax,AllKs,AllLambdas,evectors] = T_InstabilityWL(Vs,Ps,Es,varargin)
% Find the wavenumber (k) and wavelength (Lambda) of instabilities

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es.JacMode = 1;	% Request a jacobian

[vecs,evs]=eig(full(Ps.LocFunc(Vs,Ps,Es)+Ps.SpaFunc(Vs,Ps,Es)));
evsunstable=find(real(diag(evs))>0);


%[~,inds]=max(real(fft(vecs(1:Ps.Nx,evsunstable))));
%tmp=real(fft(vecs(1:Ps.Nx,evsunstable)));
%[~,inds]=max(tmp(1:min(Ps.Nx/2),:));
%K_vec = Ps.Nx/2*linspace(0,1,Ps.Nx/2+1)/Ps.Lx;
%KMinmax = K_vec(inds)*2*pi;
size(evsunstable)
for ii=1:length(evsunstable)
    %disp([ii size(vecs(:,evsunstable(ii)))]);
    tmp=T_CountRegions(reshape(real(vecs(:,evsunstable(ii))),Ps.Nx,Ps.VarNum),Ps,Es);
    AllLambdas(ii) = Ps.Lx/tmp(1);
end;
if(isempty(evsunstable))
    AllLambdas=[];
    AllKs=[];
    KMinmax=[-1 -1];
    LambdaMinmax=[-1 -1];
else
    AllKs = 2*pi./AllLambdas;
    KMinmax = [min(AllKs) max(AllKs)];
    LambdaMinmax = [min(AllLambdas) max(AllLambdas)];
end;

if(nargout>4)
    evectors=vecs(:,evsunstable);
end;

end